function[bTSONMRA,wTSONMRA,mTSONMRA,sTSONMRA,rTSONMRA,cTSONMRA]=TSO_NMRA_script
% [bFPA,wFPA,mFPA,sFPA,cFPA,rFPA]=ALL_script
n=1;%Number of runs
PopSize=10;
Iterations=100;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%FITNESS FUNCTION DETAILS%%
    Function_name='F1'
[Lb,Ub,Dim,Fun] = Get_CEC2005_Functions_details(Function_name)

%     %%%%%%%TSO-NMRA%%%%%%%%%%%%%%%%
for i=1:n
  [TSONMRAbest,TSONMRAfmin,bb]=NMRA_TSO(PopSize,Iterations,Lb,Ub,Dim,Fun);
   TSONMRAbest(i,:)=TSONMRAbest;
    rTSONMRA(i,:)=TSONMRAfmin;
    eTSONMRA(i,:)=bb;
end
disp('TSONMRA runs completed');
cTSONMRA=min(eTSONMRA);
bTSONMRA=min(rTSONMRA);
wTSONMRA=max(rTSONMRA);
mTSONMRA=mean(rTSONMRA);
sTSONMRA=std(rTSONMRA);
TSONMRAbest=min(TSONMRAbest);

%
end


