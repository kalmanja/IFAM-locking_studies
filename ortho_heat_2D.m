
    %     ________________________________
    %     |                               |<------
    %     |                               |<------
    %     |                               |<------
    %     |                               |<------
    %     |                               |<------
    %     |        2D Heat conduction     |<------Flux BC
    % T=0-|                               |<------
    % T=0-|                               |<------
    % T=0-|                               |<------
    % T=0-|                               |<------
    % T=0-|_______________________________|<------
    % T=0 0 0 0 0 0 0 0 0 0
    %           ^
    %           |
    %           |
    % Fixed temperature BC
    
    
    
    % Author: Kiran Manjunatha, IFAM (kiran.manjunatha@ifam.rwth-aachen.de)
    % Created on: 14-11-2018
    % Updated on: 16-11-2018


clear all
% clc
% clf

%% Parameters
D_1 = [1 0; 0 0.001]; % Thermal conductivity matrix
rot_ang = 20;

thk = 1;
q_bar = -1;

elem_types = 2;
mesh_sizes = 1;
corner_node_val = zeros(elem_types,mesh_sizes);
dof = zeros(elem_types,mesh_sizes);
rot_ang_rad = degtorad(rot_ang);
R = [cos(rot_ang_rad) -sin(rot_ang_rad); sin(rot_ang_rad) cos(rot_ang_rad)];
D = R*D_1*R';


for elem_type = [1 2] % 1 - four-noded isoparametric quad; 2 - eight-noded isoparametric quad 
              
    
    count = 0;
    
    if elem_type ==1
        ngauss = 4;
    elseif elem_type ==2
        ngauss = 9;
    end

    %% Mesh refinements loop
    for ne = [4 16 64 256 1024 4096]
%         length = 2;interp_length = 0.25*length/sqrt(ne);
        count = count+1;
        if elem_type == 2 && ne == 4096
            break
        end
        %% Mesh
        if elem_type == 1
            elem = importdata(strcat('elem_',num2str(ne),'.txt'));
            node = importdata(strcat('node_',num2str(ne),'.txt'));
        elseif elem_type == 2
            elem = importdata(strcat('elem_',num2str(ne),'_quadratic.txt'));
            node = importdata(strcat('node_',num2str(ne),'_quadratic.txt'));
        end
        nn = size(node,1);

        %% Stiffness matrix construction

        K = zeros(nn,nn);
        for i = 1:ne
            if elem_type==1             % 4 noded quadrilateral
                X = node(elem(i,2:5),2);
                Y = node(elem(i,2:5),3);
                nodenum = elem(i,2:5);
            elseif elem_type==2         % 8 noded quadrilateral
                X = node(elem(i,2:9),2);
                Y = node(elem(i,2:9),3);
                nodenum = elem(i,2:9);
            end
        %     X = [1;3;3.5;-1];Y = [0;-1;1;1];        % test
            kloc = stiffnessmatrix(X,Y,elem_type,ngauss,thk,D);
%             [V,e,W] = eig(kloc)
%             pause
            if elem_type == 1
                for j = 1:4
                    for k = 1:4
                        K(nodenum(j),nodenum(k)) = K(nodenum(j),nodenum(k)) + kloc(j,k);
                    end
                end
            elseif elem_type == 2
                for j = 1:8
                    for k = 1:8
                        K(nodenum(j),nodenum(k)) = K(nodenum(j),nodenum(k)) + kloc(j,k);
                    end
                end
            end
        end

        %% Boundary conditions
        [f_b,a_b] = applybc(nn,thk,q_bar);

        %% Solution
        a = zeros(nn,1);
        free_nodes = setdiff(node(:,1),a_b);
        K_F = K(free_nodes,free_nodes);
        f_b_F = f_b(free_nodes);
        a_F = K_F\f_b_F;

        a(free_nodes) = a_F;

        %% Interpolation
%         figure(ne+nn)
%         colorbar('eastoutside');
%         
%         for i = 1:ne
%             X = [node(elem(i,2),2);node(elem(i,3),2);node(elem(i,4),2);node(elem(i,5),2)];
%             Y = [node(elem(i,2),3);node(elem(i,3),3);node(elem(i,4),3);node(elem(i,5),3)];
%             aa = [a(elem(i,2));a(elem(i,3));a(elem(i,4));a(elem(i,5))];
%             patch(X,Y,aa);hold on
%         end
%         hold off

        dof(elem_type,count) = nn;
        if elem_type == 1
            corner_node_val(elem_type,count) = a(nn);
        elseif elem_type == 2
            if nn==21
                index = 9; 
            elseif nn==65
                index = 25;
            elseif nn==225
                index = 81;
            elseif nn==833
                index = 289;
            elseif nn==3201
                index = 1089;
            end
            corner_node_val(elem_type,count) = a(index);
        end
    end
    figure(1)
    if elem_type==1
        plot(dof(elem_type,:),corner_node_val(elem_type,:))
    elseif elem_type==2
        plot(dof(elem_type,1:5),corner_node_val(elem_type,1:5))
    end
    hold on
%     pause
end

