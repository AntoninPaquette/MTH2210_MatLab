function [ f ] = my_edo( t , z )

f		=	zeros(size(z)); 

f(1)	=	z(2);
f(2)	=	-10*z(1);

end

