function [ z ] = findStrcmp( x,y )

  z = find(strcmp(x,y));
  if isempty(z)
    z = 0;
  end

end

