close all; clear; clc;
%%
j       =   sqrt(-1);
%%
A       =   [ 
                1.0+j*4.0 1.0+j*2.0 1.0+j*2.0;
                1.0+j*4.0 1.0+j*3.0 1.0+j*4.0;
                1.0+j*3.0 1.0+j*4.0 1.0+j*3.0
             ];
b       =   [   1.0+j*2.0 ; 2.0+j*4.0 ; 1.0+j*4.0 ];
x       =   A\b;
A
b
x
x_      =   GaussPivot(3,A,b)
%%
function[x]=GaussPivot(N,A,b)
for k=1:N-1
    [A,b]=PivotSwap(N,A,b,k);
    if A(k,k)==0
        break;
    end
    for i=k+1:N
        b(i) = b(i)-A(i,k)*b(k)/A(k,k);
        for j=k+1:N
            A(i,j) = A(i,j)-A(i,k)*A(k,j)/A(k,k);
        end
        A(i,k)      = 0; 
    end
end
x       =   zeros(N,1);
for i=N:-1:1
    s       =   0;
    for j=i+1:N
        s       =   s+A(i,j)*x(j);
    end
    x(i,1)	=   (b(i)-s)/A(i,i);
end

%%
end
%%
function[A,b]=PivotSwap(N,A,b,k)
    Amax    =   A(k,k);
    p       =   k;
    for i=k+1:N
        if abs(A(i,k))>abs(Amax)
            Amax    =   A(i,k);
            p       =   i;
        end
    end
    if p~=k
        D       =   b(p);
        b(p)    =   b(k);
        b(k)    =   D;
        for j=1:N
            D       =   A(p,j);
            A(p,j) 	=   A(k,j);
            A(k,j) 	=   D;
        end
    end
end
%%