% y = varbirds(V)
% V, Array of vectors for velocity each row consists of coordinates for a v
function y = varbirds(V)
N = size(V,1);
onecol = ones(N,1);
Mat_dVx = V(:,1)*onecol'-onecol*(V(:,1))'; % component matrices of diffs in V
Mat_dVy = V(:,2)*onecol'-onecol*(V(:,2))';

sqrnorms = Mat_dVx.^2+Mat_dVy.^2;   % Matrix of (vj-vi)^2

y=1/(2*N^2)*sum(sum(sqrnorms));     % Calculate variance
end

%---Generalized but slow----
% norms = zeros(N);
% for i = 1:N
%     for j = 1:N
%         norms(i,j) = norm(V(i,:)-V(j,:));
%     end
% end