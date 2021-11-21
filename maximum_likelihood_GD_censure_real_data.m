 function G=maximum_likelihood_GD_censure_real_data (theta)

alpha=theta(1);
beta=theta(2);

% x=xlsread('C:\Users\Tass\Desktop\dataset Amyotrophic lateral sclerosis\pro_act_data_200_complete.xlsx');

% x=xlsread('C:\Users\Tass\Desktop\dataset Amyotrophic lateral sclerosis\pro_act_data_700_plus_censure.xlsx');
% t=x(:,1);
% censure=x(:,2);

% t=[28;8;116;2;1;112;7;17;32;1;19;5;29;42;91;8;91;85;2;1;89;3;1;5;1;84;1;112;99;
% 1;6;2;62;6;103;24;1;4;19;24;26;1;8;2;2;1;77;27;65;2;64;8;16;91;54;90;60;8;1;1;
% 93;3;101;10;5;12;1;85;37;5;89;2;91;88;44;19;1;15;69;88;25;36;86;45;7;22;90;2;
% 49;9;8;66;6;15;11;1;46;23;5;28;4;4;15;13;1;1;6;32;49;23;2;39;19;37;55;43;1;1];
% 
% censure=[0; 1; 0; 1; 0; 0; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 0; 1; 1; 1; 0; 0; 1; 0; 0; 1; 1; 1; 1; 1; 0; 1; 1;
% 0; 0; 1; 1; 0; 1; 1; 1; 1; 0; 0; 0; 1; 0; 1; 1; 0; 0; 0; 1; 1; 1; 1; 0; 0; 0; 1; 1; 1; 1; 0; 1; 1; 0; 1; 0; 0; 1; 1;
% 0; 1; 0; 0; 1; 1; 0; 0; 1; 1; 0; 1; 1; 1; 1; 0; 1; 1; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;
% 0; 0; 1];

% t=xlsread('C:\Users\Tass\Desktop\Marshall-Olkin Censure 2\real data\treatment_1.xlsx');
% censure=xlsread('C:\Users\Tass\Desktop\Marshall-Olkin Censure 2\real data\treatment_1_censure.xlsx');

% t=xlsread('C:\Users\Tass\Desktop\Marshall-Olkin Censure 2\real data\treatment_1.xlsx');
% censure=xlsread('C:\Users\Tass\Desktop\Marshall-Olkin Censure 2\real data\treatment_1_censure.xlsx');


 x=xlsread('C:\Users\Tass\Desktop\dataset Amyotrophic lateral sclerosis\pro_act_data_166_50%.xlsx');
t=x(:,1);
censure=x(:,2);


n=length(t);
 

 %% dérivé par rapport à alpha
  mc4=0;
 mc5=0;

 for i=1:n
     mc4=mc4+(censure(i));
     mc5=mc5+(exp(-beta.*t(i)));
     end
 Gompertz(1)=1./alpha.*mc4-(n-mc5)./beta;

 %% dérivé par rapport à beta
 mc10=0;
 mc11=0;
 mc12=0;
 mc13=0;
 mc14=0;
 mc15=0;
 mc16=0;
 for i=1:n
     mc10=mc10+(censure(i).*t(i));  
     mc11=mc11+((exp(-beta.*t(i))));
     mc12=mc12+(t(i).*exp(-beta.*t(i)));
       end
  Gompertz(2)=-mc10+alpha*(n-mc11)./beta^2-alpha*(mc12)/beta;
%   Gompertz(3)=beta;

 %%
G=[Gompertz(1), Gompertz(2)];
 end
