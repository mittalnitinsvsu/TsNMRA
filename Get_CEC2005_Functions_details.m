
% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim,fobj] = Get_CEC2005_Functions_details(F)


switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=200;
        
    case 'F2'
        fobj = @F2;
        lb=-10;
        ub=10;
        dim=200;
        
              
end

end

% F1

function o = F1(x)
o=sum(x.^2);
end

% F2

function o = F2(x)
o=sum(abs(x))+prod(abs(x));
end


function o=Ufun(x,a,k,m)
o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end