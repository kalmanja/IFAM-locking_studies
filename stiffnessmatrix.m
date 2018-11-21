function kloc = stiffnessmatrix(X,Y,elem_type,ngauss,thk,D)

if ngauss == 1
    [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,0,0);
    kloc = B'*thk*D*B*jac;
end

if ngauss == 4
    kloc = zeros(4,4);
    for i = 1:ngauss
        if i == 1
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,-1/sqrt(3),-1/sqrt(3));
        elseif i ==2
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,1/sqrt(3),-1/sqrt(3));
        elseif i == 3
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,1/sqrt(3),1/sqrt(3));
        else
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,-1/sqrt(3),1/sqrt(3));
        end
        kloc = kloc + B'*thk*D*B*jac;
    end
end

if ngauss == 9
    kloc = zeros(8,8);
    for i = 1:ngauss
        if i == 1
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,-sqrt(3/5),-sqrt(3/5));
            kloc = kloc + (25/81)*B'*thk*D*B*jac;
        elseif i ==2
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,0,-sqrt(3/5));
            kloc = kloc + (40/81)*B'*thk*D*B*jac;
        elseif i == 3
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,sqrt(3/5),-sqrt(3/5));
            kloc = kloc + (25/81)*B'*thk*D*B*jac;
        elseif i == 4
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,-sqrt(3/5),0);
            kloc = kloc + (40/81)*B'*thk*D*B*jac;
        elseif i == 5
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,0,0);
            kloc = kloc + (64/81)*B'*thk*D*B*jac;  
        elseif i == 6
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,sqrt(3/5),0);
            kloc = kloc + (40/81)*B'*thk*D*B*jac;
        elseif i == 7
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,-sqrt(3/5),sqrt(3/5));
            kloc = kloc + (25/81)*B'*thk*D*B*jac;
        elseif i == 8
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,0,sqrt(3/5));
            kloc = kloc + (40/81)*B'*thk*D*B*jac;
        elseif i == 9
            [J,Jinv,jac,B] = jacobian_calc(elem_type,X,Y,sqrt(3/5),sqrt(3/5));
            kloc = kloc + (25/81)*B'*thk*D*B*jac;
        end
    end
end
end