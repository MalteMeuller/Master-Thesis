%converts vector to symmetrix matrix
function VV=unvech(vv)

m1=((1+8*length(vv))^0.5-1)/2;


VV=zeros(m1,m1);
counter1=1;
%populates to matrix
for im1=1:m1
    for im2=im1:m1
        VV(im2,im1)=vv(counter1,1);
        counter1=counter1+1;
    end
end

for im1=1:m1
    for im2=im1:m1
        VV(im1,im2)=VV(im2,im1);
    end
end
