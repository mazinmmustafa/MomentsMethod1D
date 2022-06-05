function[y]=sinc(x)
%%
y       =	sin(x)./(x);
y(~x)	=	1;
end
%%