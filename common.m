# Prevent Octave from thinking that this is a script file to execute
1;

function draw_line(x1, y1, x2, y2, c = "black")
# Usage: draw_line(x1, y1, x2, y2, c = "black")
#
# Draws a line from the coordinate x1/y1 to x2/y2
# with optional color.

    line([x1,x2],[y1,y2], "color", c);
endfunction

function draw_horizontal_line(x1, x2, y, c = "black")
# Usage: draw_horizontal_line(x1, x2, y, c = "black")
#
# Draws a horizontal line at y at coordinates from x1 to x2
# with optional color.

    line([x1,x2],[y,y], "color", c);
endfunction

function draw_vertical_line(x, y1, y2, c = "black")
# Usage: draw_vertical_line(x, y1, y2, c = "black")
#
# Draws a vertical line at coordinate x from y1 to y2
# with optional color.

    line([x,x],[y1,y2], "color", c);
endfunction

function draw_cross_line(x, y, c = "red", scale = 40)
# Usage: draw_cross_line(x, y, c = "red", scale = 40)
#
# Draws a cross line at coordinate x/y
# with optional color and cross size scaling.
    xl = get (gca(), "xlim");
    h = (xl(2)-xl(1)) / scale;
    yl = get (gca(), "ylim");
    v = (yl(2)-yl(1)) / scale;
    draw_horizontal_line(x-h, x+h, y, c);
    draw_vertical_line(x, y-v, y+v, c);
endfunction    

function draw_circle(xc, yc, r, c = "blue", step=100)
# Usage: draw_cirle(xc, yc, r, step=100, c = "red")
#
# Draws a circle at center xc/yc with radius r and optional color and step size
    theta = linspace(0, 2*pi, step);
    x = r * cos(theta) + xc;
    y = r * sin(theta) + yc;
    plot(x, y, "color", c, "Linewidth", 2);
endfunction   
    
function y = fill_vector(cnt, val = 0)
# Usage: y = fill_vector(cnt, val = 0)
#
# Returns a vector 'y' with 'cnt' elements of value 'val'.

    y = val(ones(1,cnt));
endfunction

function idx = idx_min(y)
# Usage: idx = idx_min(y)
#
# Returns the index 'idx' within 'y' holding the minimum value.

    [v, i] = min(y);
    idx = i;
endfunction

function idx = idx_max(y)
# Usage: idx = idx_max(y)
#
# Returns the index 'idx' within 'y' holding the maximum value.

    [v, i] = max(y);
    idx = i;
endfunction

function d = diff2(y, r = 0)
# Usage: d = diff2(y, r = 0)
#
# Computes the differiantial equation of y with radius r.

    if (r == 0)
        d = [0, diff(y)];
        return;
    endif

    [row, col] = size(y);
    dd = zeros(1,col-2*r);
    dd = y(1+2*r:end) - y(1:end-2*r);
    d = [zeros(1,r), dd, zeros(1,r)];
endfunction

function [par, fmin, info] = least_square(f_fit, x, y, p0)
# Usage: [par, fmin, info] = least_square(f_fit, x, y, p0)
#
# Solves function approximation with least square approach.

# Least square loss function
f_ls = @(p) sum((f_fit(x,p)-y).^2);

# Parameter optimization
[par, fmin, info] = fminunc(f_ls, p0);
endfunction

function l = mse(y)
# Usage: l = mse(y)
#
# Returns a vector 'l' with the sliding mean square error from
# vector 'y'.

    s = numel(y);
    l = fill_vector(s);
    
    for i = 2:s-1
        for j = 2:i 
            l(i) += ( y(j) - mean(y(1:j)) )^2; 
         endfor; 
         for j = i+1:s-1
            l(i) += ( y(j) - mean(y(j:s)) )^2;
        endfor
    endfor
    
    # Fill ignored borders with their neighbor values
    l(1) = l(2);
    l(s) = l(s-1);
endfunction

function sm = submatrix(M, r, c, radius)
# Usage: sm = submatrix(M, r, c, radius)
#
# Returns a submatrix 'sm' out of bigger matrix surround center row 'r'/column 'c'
# within "radius".

    sm = M(-radius+r:r+radius, -radius+c:c+radius);
endfunction

function fn_list = get_files_from_dir(d, fp)
# Usage: fn_list = get_files_from_dir(d, fp)
#
# Returns a file name list from a given directory d according file pattern fp.

    fn_list = glob(strcat(d, fp));
endfunction

function convert_tro_csv(fn_in)
# Usage: convert_tro_csv(fn_in)
#
# Convert output of TR offline csv in standard profile csv and save it.

    s = strsplit(fn_in, ".");
    fn_out = strcat(s{1}, "-ptb.", s{2});

    M = dlmread(fn_in, ";");

    r = columns(M);
    printf("Found matrix with columns %i\n", r);
 
    N(1,:) = [0:r-1];
    
    #N(2,:) = M(2, :) * 10;   # Height 10 um to um // Old version of TRO
    #N(3,:) = M(3, :) * -10;  # Height 10 um to um
    
    N(2,:) =  M(2, :);  # Height in um
    N(3,:) = -M(3, :);  # Height in um
 
    plot(N(1,:), N(2,:), N(3,:));
    grid on;
    
    dlmwrite (fn_out, N, ";");
endfunction

function convert_bf_csv(fn_in, length = 200, mirror = false)
# Usage: convert_bf_csv(fn_in)
#
# Convert output of beadfiller csv in standard profile csv and save it.

    s = strsplit(fn_in, ".");
    fn_out = strcat(s{1}, "-ptb.", s{2});
    
    m = dlmread(fn_in, ";");
    
    if (rows(m) > 1)
        error("Only vector as input is allowed!");
    endif
    
    turn_idx = 0;
    i=1;
    while (i <= columns(m)-2)
        if(m(i+2) < m(i))
            turn_idx = i;
            printf("Found turn point at i/m(i)/m(i+1) %i/%f/%f\n", turn_idx, m(turn_idx), m(turn_idx+1));
            break;
        endif;
        i = i+2;
    endwhile;
    
    m1 = m(1:turn_idx-1);
    m2 = m(turn_idx:end);
    
    printf("Found vector length m1/m2 %i/%i\n", columns(m1), columns(m2));
        
    x1 = round(m1(1:2:end) * 10);   # Position to index
    z1 = round(m1(2:2:end) * 1000); # Height mm to um
    
    printf("Vector length x1/z1 %i/%i\n", columns(x1), columns(z1));
        
    if (columns(x1) != columns(z1))
        error("Vector size of x1 and z1 do not match!");
    endif
    
    x2 = round(flip(m2(1:2:end)) * 10);   # Position to index
    z2 = round(flip(m2(2:2:end)) * 1000); # Height mm in um
    
    printf("Vector length x2/z2 %i/%i\n", columns(x2), columns(z2));
    
    if (columns(x2) != columns(z2))
        error("Vector size of x2 and z2 do not match!");
    endif
        
    if (columns(x1) >= columns(x2))
        xx1 = x1(1:columns(x2));
        zz1 = z1(1:columns(x2));
        xx2 = x2;
        zz2 = z2;
    else
        xx1 = x1;
        zz1 = z1;
        xx2 = x2(1:columns(x1));
        zz2 = z2(1:columns(x1));
    endif
    
    if (mirror)
        zz1 = flip(zz1);
        zz2 = flip(zz2);
    endif
    
    printf("Found point set at x1/z1 (%i/%i) (%i/%i) and x2/z2 (%i/%i) (%i/%i)\n",
           xx1(1), zz1(1), xx1(end), zz1(end), xx2(1), zz2(1), xx2(end), zz2(end));
    
    subplot (1, 2, 1);
    plot(xx1, zz1, xx2, zz2);
    grid on;
    ylabel ("height [um]");
    
    M(1,:) = xx1(1:length);
    M(2,:) = zz1(1:length);
    M(3,:) = zz2(1:length);
    
    subplot (1, 2, 2);
    plot(M(1,:), M(2,:), M(1,:), M(3,:));
    grid on;
    ylabel ("height [um]");
    
    dlmwrite (fn_out, M, ";");
endfunction

function convert_sol_csv(fn_in)
# Usage: convert_sol_csv(fn_in)
#
# Convert output of SOL csv in standard profile csv and save it.

    s = strsplit(fn_in, ".");
    M = dlmread(fn_in, ",");
    
    r = rows(M);
    printf("Found matrix with rows %i\n", r);
    c = columns(M);
    printf("Found matrix with columns %i\n", c);
    
    fig_no = 0;
    
    for i = 1:r
        fn_out = strcat(s{1}, "-", num2str(i), "-ptb.", s{2});

        N(1,:) = [0:c-1];
        N(2,:) = M(i, :);       # Top, height in um
        N(3,:) = zeros(1,c);    # Bottom, set to zero
        
        figure(++fig_no);
        semilogy(N(1,:), N(2,:));
        grid on;
    
        dlmwrite (fn_out, N, ";");
    endfor
endfunction

function p = laser_tracking_cog(I, radius)
# Usage: p = laser_tracking_cog(I, radius)
#
# Returns a vector 'p' with the laser line column indices from
# matrix 'I' as intensive image. Starts laser line tracking from
# max. intensive in the column within 'radius'.
    
    r = rows(I);
    c = columns(I);
    ri = zeros(1,c);
    p = zeros(1,c);

    for j = [1:c]
        [v, vi] = max(I(:, j));
        if vi <= radius
            ri(j) = 1+radius;
            continue;
        elseif vi > r-radius;
            ri(j) = r-radius-1;
            continue;
        endif
        ri(j) = vi;
    endfor
    
    for j = [1:c]
        iv = [-radius+ri(j):ri(j)+radius];
        vv = double(I(-radius+ri(j):ri(j)+radius, j));
        s = sum(vv);
        if s > 1
            p(j) = (iv*vv) / s;
        else
            p(j) = 1;
        endif
    endfor
endfunction

function [value8, length] = read_binary(filename)
# Usage: [value8, length] = read_binary(filename)
#
# Returns a raw vector 'value8' of type uint8 and its length
# from content of file 'filename'.
	fid = fopen(filename);
	[data, length] = fread(fid);
    value8 = transpose(uint8(data));
	fclose(fid);
endfunction

function [length, profile, footer] = convert_jsonplus(value8)
# Usage: [length, profile, footer] = convert_jsonplus(value8)
#
# Returns the 'length' of raw vector 'profile' of type int32 [um] and
# the 'footer' as JSON string.
	length = 0;
    profile = [];
    footer = [];
    
    # Extract the lenght of binary part
    _10_7 = double(value8(1)) - 0x30;
    _10_6 = double(value8(2)) - 0x30;
    _10_5 = double(value8(3)) - 0x30;
    _10_4 = double(value8(4)) - 0x30;
    _10_3 = double(value8(5)) - 0x30;
    _10_2 = double(value8(6)) - 0x30;
    _10_1 = double(value8(7)) - 0x30;
    _10_0 = double(value8(8)) - 0x30;
    
    length = _10_7 * 1E7 + _10_6 * 1E6 + _10_5 * 1E5 + _10_4 * 1E4 +...
             _10_3 * 1E3 + _10_2 * 1E2 + _10_1 * 1E1 + _10_0 * 1E0;

    # Extract profile
    for i = [1:length/4]
        j = 5 + 4 * i;
        profile(i) = typecast([value8(j), value8(j+1), value8(j+2), value8(j+3)], 'int32');
    endfor    
    
    # Extract JSON footer
    footer = char(value8(9+length:end));
endfunction





