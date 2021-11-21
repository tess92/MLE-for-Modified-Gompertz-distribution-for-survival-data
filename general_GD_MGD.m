clear all
close all

alpha=0.003; 
beta=0.0001; 


%%% il faut que beta soit negative pour qu'il soit defective (spécialement pour ce programme)


options=optimset('disp','iter','LargeScale','off','TolFun',.001,'MaxIter',10000,'MaxFunEvals',10000);

fsolve(@maximum_likelihood_GD_censure_real_data , [ alpha beta]')