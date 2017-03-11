function fitsed_fitfunc, x, m

basis = dblarr(m)

for i=0,m-1 do basis[i] = x

return, basis

end

