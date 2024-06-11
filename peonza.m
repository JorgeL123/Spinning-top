function [DDT,DDF,DDP]=peonza(t,f,p,dt,df,dp,I1,I3,h,M);
g=9.81;
DDT=(df^2*sin(t)*cos(t)*(I1-I3)-I3*(dp)*df*sin(t)+M*g*h*sin(t))/I1;
DDF=(2*(I3-I1)*df*cos(t)*dt*sin(t)-I3*(df*cos(t))*sin(t)*dt+I3*dp*dt*sin(t))/(I1*(sin(t))^2);
DDP=(cos(t)*((I1*sin(t)^2+I3*cos(t)^2)*df*dt*sin(t)/cos(t)-2*(I3-I1)*df*dt*sin(t)*cos(t)-I3*dp*dt*sin(t))/(I1*(sin(t))^2));
end