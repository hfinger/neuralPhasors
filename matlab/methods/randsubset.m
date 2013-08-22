function [r] = randsubset(n,k)

assert(n>=k)

if k>n/2
  k=n-k;
  inverted=true;
else
  inverted=false;
end

r = randi(n,[1 k]);
r = unique(r);
while length(r)<k
  r = union(r, randi(n,[1 k-length(r)]));
end

if inverted
  r = setdiff(1:n,r);
end

end