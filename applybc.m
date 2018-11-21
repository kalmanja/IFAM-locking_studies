function [f_b,a_b] = applybc(nn,thk,q_bar) 

f_b = zeros(nn,1);
        if nn==9
            % 4 elements - linear
            % flux BC
            gridlength = 1.0;
            f_b(3) = -gridlength*thk*q_bar/2;
            f_b(6) = 2*f_b(3);
            f_b(9) = f_b(3);
            a_b = [1 2 4];  % Essential BC (T = 0) at defined nodes
        elseif nn==21
            % 4 elements - quadratic
            % flux BC
            gridlength = 1.0;
            f_b(3) = -gridlength*thk*q_bar/6;
            f_b(6) = 2*f_b(3);
            f_b(9) = f_b(3);
            f_b([15 20])=4*f_b(3);
            a_b = [1 2 4 10 13];  % Essential BC (T = 0) at defined nodes
        elseif nn==25
            % 16 elements  - linear
            gridlength = 0.5;
            f_b([5  25]) = -gridlength*thk*q_bar/2;
            f_b([10 15 20]) = -gridlength*thk*q_bar;
            a_b = [1    2   3   6   11];  % Essential BC (T = 0) at defined nodes
        elseif nn==65
            % 16 elements - quadratic
            % flux BC
            gridlength = 0.5;
            f_b(5) = -gridlength*thk*q_bar/6;
            f_b([10 15 20]) = 2*f_b(5);
            f_b(25) = f_b(5);
            f_b([37 46 55 64])=4*f_b(5);
            a_b = [1 2 3 6 11 26 29 30 41];  % Essential BC (T = 0) at defined nodes
        elseif nn==81
            % 64 elements - linear
            gridlength = 0.25;
            f_b([9  81]) = -gridlength*thk*q_bar/2;
            f_b(18:9:72) = -gridlength*thk*q_bar;
            a_b = [1    2   3   4   5   10  19  28  37];  % Essential BC (T = 0) at defined nodes
        elseif nn==225
            % 64 elements - quadratic
            % flux BC
            gridlength = 0.25;
            f_b(9) = -gridlength*thk*q_bar/6;
            f_b(18:9:72) = 2*f_b(9);
            f_b(81) = f_b(9);
            f_b(105:17:224)=4*f_b(9);
            a_b = [1:5 10:9:37 82 86:3:92 85 109 126 143];  % Essential BC (T = 0) at defined nodes
        elseif nn==289
            % 256 elements - linear
            gridlength = 0.125;
            f_b([17  289]) = -gridlength*thk*q_bar/2;
            f_b(34:17:272) = -gridlength*thk*q_bar;
            a_b = [1    2   3   4   5   6   7   8   9   18  35  52  69  86  103 120 137];  % Essential BC (T = 0) at defined nodes
        elseif nn==833
            % 256 elements - quadratic
            % flux BC
            gridlength = 0.125;
            f_b(17) = -gridlength*thk*q_bar/6;
            f_b(34:17:272) = 2*f_b(17);
            f_b(289) = f_b(17);
            f_b(337:33:832)=4*f_b(17);
            a_b = [1:9 18:17:137 290 294:3:312 293 341:33:539];  % Essential BC (T = 0) at defined nodes
        elseif nn==1089
            % 1024 elements - linear
            gridlength = 0.0625;
            f_b([33  1089]) = -gridlength*thk*q_bar/2;
            f_b(66:33:1056) = -gridlength*thk*q_bar;
            a_b = [1:17   34:33:529];  % Essential BC (T = 0) at defined nodes
        elseif nn==3201
            % 1024 elements - quadratic
            % flux BC
            gridlength = 0.0625;
            f_b(33) = -gridlength*thk*q_bar/6;
            f_b(66:33:1056) = 2*f_b(33);
            f_b(1089) = f_b(33);
            f_b(1185:65:3200)=4*f_b(33);
            a_b = [1:17 34:33:529 1090 1094:3:1136 1093 1189:65:2099];  % Essential BC (T = 0) at defined nodes
        elseif nn==4225
            % 1024 elements - linear
            gridlength = 0.03125;
            f_b([65  4225]) = -gridlength*thk*q_bar/2;
            f_b(130:65:4160) = -gridlength*thk*q_bar;
            a_b = [1:33   66:65:2081];  % Essential BC (T = 0) at defined nodes
        end
        
end