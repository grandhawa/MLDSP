% Voronoi polygons:
load Griffith_Peres-Neto_2006.mat

x = xy(:,1);
y = xy(:,2);
z = zeros(1,size(x,1));

tri = delaunay(x,y);

hold on;
trimesh(tri,x,y,z);
% hidden off;


% plot(x,y,'k.');

voronoi(x,y,tri);

[v,c] = voronoin([x y]);

% plot(v(2:end,1),v(2:end,2),'r.')

for i = 1:size(c,1)

   plot(v(c{i},1),v(c{i},2),'-r.')

end



k = dsearch(x,y,tri,x(3),y(3))