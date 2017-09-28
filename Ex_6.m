%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of Examples from TOOLS FOR THE STUDY OF STABILITY AND CONVERGENCE IN SET
% DYNAMICAL SYSTEMS WITH APPLICATIONS TO FEEDBACK CONTROL
% Example: 6
% Nathalie Risso. nrisso@email.arizona.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1000; % set density
X0=linspace(0.5,1,n);
J=10;
B=1;
X=zeros(n,J);
A = [1]; 
dmax = 0.2;
ball = linspace(-1,1,n);
X(:,1)=X0;
k = -0.6;
F = A +B*k+dmax*ball;
F=F';
for j=1:J-1
    temp = allcomb((F),(X(:,j)));
    tp = resample(temp,1,1000);
     X(:,j+1)=tp(:,1).*tp(:,2);
   end
figure(1)
hold on;
for i=1:J
   plot(0*X(:,i)+i-1,X(:,i),'b','LineWidth',2); 
  xlabel('j', 'Interpreter','latex')
   hold on;
   grid on
end
