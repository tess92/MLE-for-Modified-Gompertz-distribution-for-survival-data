 clear all
close all
clc


% alpha=0:0.0001:1;
% beta=0:0.0001:0.001; 

alpha=0:0.0001:0.01;
beta=-0.01:0.0001:0;

x=xlsread('C:\Users\Tass\Desktop\dataset Amyotrophic lateral sclerosis\pro_act_data_700_plus_censure.xlsx');
t=x(:,1);
censure=x(:,2);


% t=[28;8;116;2;1;112;7;17;32;1;19;5;29;42;91;8;91;85;2;1;89;3;1;5;1;84;1;112;99;
% 1;6;2;62;6;103;24;1;4;19;24;26;1;8;2;2;1;77;27;65;2;64;8;16;91;54;90;60;8;1;1;
% 93;3;101;10;5;12;1;85;37;5;89;2;91;88;44;19;1;15;69;88;25;36;86;45;7;22;90;2;
% 49;9;8;66;6;15;11;1;46;23;5;28;4;4;15;13;1;1;6;32;49;23;2;39;19;37;55;43;1;1];
% censure=[0; 1; 0; 1; 0; 0; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 0; 1; 1; 1; 0; 0; 1; 0; 0; 1; 1; 1; 1; 1; 0; 1; 1;
% 0; 0; 1; 1; 0; 1; 1; 1; 1; 0; 0; 0; 1; 0; 1; 1; 0; 0; 0; 1; 1; 1; 1; 0; 0; 0; 1; 1; 1; 1; 0; 1; 1; 0; 1; 0; 0; 1; 1;
% 0; 1; 0; 0; 1; 1; 0; 0; 1; 1; 0; 1; 1; 1; 1; 0; 1; 1; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;
% 0; 0; 1];


% t=xlsread('C:\Users\Tass\Desktop\Marshall-Olkin Censure 2\real data\treatment_1.xlsx');
% censure=xlsread('C:\Users\Tass\Desktop\Marshall-Olkin Censure 2\real data\treatment_1_censure.xlsx');

n=length(t);

[xx,yy]=meshgrid(alpha,beta);





mc4=0; mc10=0; mc11=0; mc18=0; mcc=0; mccc=0; mcccc=0;
for i=1:n
mc4=mc4+(censure(i));
 mc10=mc10+(censure(i).*t(i));  
mc11=mc11+(exp(-yy.*t(i)));

end
ll=log(xx).*mc4-yy*mc10-xx.*(n-mc11)./yy;
max(ll(:))
figure
% C = gradient(ll);
mesh(xx,yy,ll)
% plot(alpha,l)
xlabel('alpha')
ylabel('beta')
zlabel('log-likelihood function')
% shading interp
