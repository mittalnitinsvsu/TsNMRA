function a=myarray(Sol,i,r)
% myarray=[1,2,3,4,5; 2,3,4,5,6; 3,4,5,6,7; 4,5,6,7,8; 5,6,7,8,9; 6,7,8,9,10; 7,8,9,10,11];
b = size(Sol);
col = b(2);
row = b(1);
% i=6;
n=r;
a = zeros(2*n+1,1);
j = 1;
while j<=2*n+1
    a(j)=i-n+j-1;
    if a(j)<1
        a(j) = a(j)+row;
    end
    if a(j)>row
        a(j) = a(j)-row;
    end
    j=j+1;
end
j=1;
while j<=2*n+1
    a;
    j=j+1;
end