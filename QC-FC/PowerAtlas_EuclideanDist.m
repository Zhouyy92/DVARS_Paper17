


for nx=1:264
    for ny=1:264
        c0=Power2011vx(nx,:);
        c1=Power2011vx(ny,:);
        ED(nx,ny)=sqrt(sum((c0-c1).^2));
    end
end