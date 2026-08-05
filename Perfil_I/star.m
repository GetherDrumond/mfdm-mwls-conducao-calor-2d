function [S] = star(N,x,y,X,Y,m)
S = [1:N; sqrt((x-X').^2+(y-Y').^2)]';
S = sortrows(S,2);
S = S(1:m,:);
