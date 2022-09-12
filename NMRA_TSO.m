%___________________________________________________________________%
% Transient Search Naked Mole-rat optimizer(NMRA) source codes       %
%            By "Nitin Mittal"                   %
%___________________________________________________________________%

function [NMRbest,fmin,bb]=NMRA_TSO(n,maxiter,Lb,Ub,d,Fun)
display('NMRATSO Working');
counter=0;
bp=0.05;
breeders=n/5;        % breeder population
% workers=n-breeders;  % worker population
popmin=1;
popmax=n;
iter=1;

for i=1:n,
    %     Lb+(Ub-Lb).*rand(1,d)
    NMRsolution(i,:)=Lb+(Ub-Lb).*rand(1,d);
    NMRfitness(i)=Fun(NMRsolution(i,:));

end

[fmin,I]=min(NMRfitness);
NMRbest=NMRsolution(I,:);
S=NMRsolution;
% W=zeros(n,d);

while (iter<=maxiter) % Loop over worker and breeders
    %For workers
        t=2-iter*(2/maxiter); 
        K=1;   % K is a real can be 0, 1, 2,....
                r1=rand();r2=rand(); 
        r3 = rand();
        T=2*t*r1-t;  
        C1=K*r2*t+1;
        alpha_min=0.5;%0.25;
        alpha_max=0.9; k=rand; p=0.95;
        lambda=alpha_min+((alpha_max-alpha_min)*p^(k-1));
%         lambda=rand;
    for i=((n/5)+1):n,
        
        %         lambda=2-iter*((2)/maxiter);
        % Find random NMR in the neighbourhood
        r=floor(randi((popmax/2)-popmin));
                 
                 
                 
                 ind=myarray(S,i,r);
                 
                 
                 NMRfitnessLNS=[NMRfitness(ind)' ind];
                 [aa,bb]=min(NMRfitnessLNS(:,1));%%min
                 
                 S_opt=NMRfitnessLNS(bb);
                 p_r=min(ind);
                 q_r=max(ind);
                 abc=randi(q_r-p_r);
                 def=randi(q_r-p_r);
                 mm=rand;
                 nn=rand;        
                 S(i,:)=S(i,:)+mm*(S(ind(bb),:)-S(i,:))+nn*(S(abc+p_r,:)-S(def+p_r,:)); 
        ab=randperm(n);
        L=Levy(d);
        %         if bp <0.05,
        if iter<maxiter/2
            S(i,:)=(NMRsolution(i,:)+L.*(NMRsolution(ab(1),:)-NMRsolution(ab(2),:)));
        else
               if r3<0.5
           
               S(i,:)=NMRbest+exp(-T)*(S(i,:)-C1*NMRbest);
                
               elseif r3>=0.5
            
               S(i,:)=NMRbest+exp(-T)*(cos(T*2*pi)+sin(T*2*pi))*abs(S(i,:)-C1*NMRbest);

              end  
        end
    end
             
        Fnew=Fun(S(i,:));
        % If NMRfitness improves (better NMRsolutions found), update then
        if (Fnew<=NMRfitness(i)),
            NMRsolution(i,:)=S(i,:);
            NMRfitness(i)=Fnew;
        end

    %For Breeders
    for z=1:breeders;
        if rand>bp
            %                 stepsize(i,j)=RB(i,j)*(Elite(i,j)-RB(i,j)*Prey(i,j));
            %                 Prey(i,j)=Prey(i,j)+P*R*stepsize(i,j);
            NMRneighbours=randperm(breeders);
            S(z,:)=(1-lambda).*(S(z,:))+(lambda*(NMRbest-NMRsolution(NMRneighbours(1),:)));
%                         Fnew= feval(Fun,S(z,:)',func_num)';
%             S(z,:)=(1/4)*sin(pi*(S(z,:)));
%                 S(z,:)=cos(0.5*acos(S(z,:)));
            Fnew=Fun(S(z,:));
            % If NMRfitness improves (better NMRsolutions found), update then
            if (Fnew<=NMRfitness(z)),
                NMRsolution(z,:)=S(z,:);
                NMRfitness(z)=Fnew;
            end
        end
    end
    if  Fnew == NMRfitness,
        counter=counter+1;
    end
    for i=1:n
        Flag4Ub=S(i,:)>Ub';
        Flag4Lb=S(i,:)<Lb';
        S(i,:)=(S(i,:).*(~(Flag4Ub+Flag4Lb)))+Ub'.*Flag4Ub+Lb'.*Flag4Lb;
    end
    %     NMRfitness=sort(NMRfitness);
    [fmin,I]=min(NMRfitness);
    NMRbest=NMRsolution(I,:);
    S=NMRsolution;
    %     S=(1/4)*sin(pi*S);
    bb(iter)=fmin;
    iter=iter+1;
    if counter <=10;
        beta=1;
        sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
        for j=1:n,
            ab=randperm(n);
            a=2-j*((2)/maxiter); % a decreases linearly fron 2 to 0
            s=S(j,:);
            u=randn(size(s))*sigma;
            v=randn(size(s));
            step=u./abs(v).^(1/beta);
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2;% Equation (3.4)
            D_alpha=abs(C1*NMRbest-S(j,:)); % Equation (3.5)-part 1
            X1=NMRbest-A1*D_alpha;% Equation (3.6)-part 1
            r1=rand();
            r2=rand();
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            D_beta=abs(C2*NMRbest-S(j,:)); % Equation (3.5)-part 2
            X2=NMRbest-A2*D_beta; % Equation (3.6)-part 2
            r1=rand();
            r2=rand();
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            D_delta=abs(C3*NMRbest-S(j,:)); % Equation (3.5)-part 3
            X3=NMRbest-A3*D_delta; % Equation (3.5)-part 3
            S(j,:)=(X1+X2+X3)/3;% Equation (3.7)
            %  s=S(j,:);
            %  stepsize=0.01*step.*(s-NMRbest);
            %  S(j,:)=s+stepsize.*randn(size(s));
            %      S(j,:)(j,:)=s;
            Fnew=Fun(S(j,:));
            %     % If NMRfitness improves (better NMRsolutions found), update then
            if (Fnew<=NMRfitness(j)),
                NMRsolution(z,:)=S(z,:);
                NMRfitness(z)=Fnew;
            end
        end
        for i=1:n
            Flag4Ub=S(i,:)>Ub';
            Flag4Lb=S(i,:)<Lb';
            S(i,:)=(S(i,:).*(~(Flag4Ub+Flag4Lb)))+Ub'.*Flag4Ub+Lb'.*Flag4Lb;
        end
    end
end

function L=Levy(d)
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
L=0.01*step;
