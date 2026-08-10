function [ jac ] = my_sys_nl_jac( x )

jac(1,1)	=	2*x(1);
jac(1,2)	=	2*x(2);
jac(2,1)	=	-2*x(1);
jac(2,2)	=	1;

end