function a_interp = postprocess(a,length,interp_length,elem,node)
map = [1 -1 -1 1; 1 1 -1 -1; 1 1 1 1; 1 -1 1 -1];
xx = 0:interp_length:length;
a_interp = zeros(size(xx,2));
[XX,YY] = meshgrid(xx,xx);
nn = size(node,1);
ne = size(elem,1);
dist = zeros(nn,1);
for i = 1:size(xx,2)
    for j = 1:size(xx,2)
       for k = 1:nn
           dist(k) = sqrt((XX(i,j)-node(k,2))^2 + (YY(i,j)-node(k,3))^2);
       end
       [sorted_dist,indices] = sort(dist,'ascend');
       nearest_nodes = indices(1:4);
       sorted_nearest_nodes = sort(nearest_nodes,'ascend');
       if sorted_dist(1) == 0
           a_interp(i,j) = a(indices(1));
       else           
           for l = 1:ne
               nodes = elem(l,2:5)';
               sorted_nodes = sort(nodes,'ascend');
               found = isequal(sorted_nearest_nodes,sorted_nodes);
               if found == 1
                    X = [node(elem(l,2),2);node(elem(l,3),2);node(elem(l,4),2);node(elem(l,5),2)];
                    Y = [node(elem(l,2),3);node(elem(l,3),3);node(elem(l,4),3);node(elem(l,5),3)];
                    invmapX = map\X;
                    invmapY = map\Y;
                    if invmapX(2)==0 || invmapX(4)==0
                        eta = (invmapX(1)-XX(i,j))/invmapX(3);
                        xi = ((invmapY(1)-YY(i,j))*invmapX(3)-(invmapX(1)-XX(i,j))*invmapY(3))/(invmapX(3)*invmapY(2) + (invmapX(1)-XX(i,j))*invmapY(4));
                        N = shapefunctionmatrix(xi,eta);
                        a_loc = [a(elem(l,2));a(elem(l,3));a(elem(l,4));a(elem(l,4))];
                        a_interp(i,j) = N*a_loc;
                        break
                    else
                        eta = -invmapX(2)/invmapX(4);
                        xi = ((invmapY(1)-YY(i,j))*invmapX(4)+invmapX(2)*invmapY(3))/(invmapX(4)*invmapY(2) - invmapX(2)*invmapY(4));
                        N = shapefunctionmatrix(xi,eta);
                        a_loc = [a(elem(l,2));a(elem(l,3));a(elem(l,4));a(elem(l,4))];
                        a_interp(i,j) = N*a_loc;
                        break
                    end
               end
           end
       end
    end
end

end
