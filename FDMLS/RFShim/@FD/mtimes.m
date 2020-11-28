function res = mtimes(a,b)

if a.adjoint
    res = adjDx(b(:,:,:,1)) + adjDy(b(:,:,:,2));
else
    res = D(b);
end

function y = D(x)

Dx = x([1,1:(end-1)],:,:) - x([2:end,end],:,:);
Dy = x(:,[1,1:(end-1)],:) - x(:,[2:end,end],:);

y = cat(4,Dx,Dy);

function y = adjDx(x)
y = x([1,1:end-1],:,:) - x;
y(1,:,:) = -x(1,:,:);
y(end,:,:) = x(end-1,:,:);

function y = adjDy(x)
y = x(:,[1,1:end-1],:) - x;
y(:,1,:) = -x(:,1,:);
y(:,end,:) = x(:,end-1,:);
