clc;
clear;
close all;
dt=0.005;
tf=9;
t=[0:dt:tf];
h=0.05;Ra=0.03;
M=0.005;I1=M*(3*h^2/5+3*Ra^2/20); I3=3*M*Ra^2/10; hcm=3*h/4;
g=9.81;
df(1)=1.5*pi;dth(1)=0;th(1)=pi/4;f(1)=0;p(1)=0; dp(1)=132*pi; %10, 120pi / 40, 90pi/ 6, 120pi/ 0, 120pi/ 1.5pi, 132pi
[ddth(1),ddf(1),ddp(1)]=peonza(th,f,p,dth,df,dp,I1,I3,hcm,M);
N=length(t);
for i=1:N-1;
    kth1=dth(i);
    kf1=df(i);
    kp1=dp(i);
        dkth1=ddth(i);
    dkf1=ddf(i);
    dkp1=ddp(i);
    %
    kth2=kth1+dkth1*dt/2;
    kf2=kf1+dkf1*dt/2;
    kp2=kp1+dkp1*dt/2;
    [dkth2,dkf2,dkp2]=peonza(th(i)+kth1*dt/2,f(i)+kf1*dt/2,p(i)+kp1*dt/2,dth(i)+dkth1*dt/2,df(i)+dkf1*dt/2,dp(i)+dkp1*dt/2,I1,I3,hcm,M);
    %
   kth3=kth1+dkth2*dt/2;
    kf3=kf1+dkf2*dt/2;
    kp3=kp1+dkp2*dt/2;
    [dkth3,dkf3,dkp3]=peonza(th(i)+kth2*dt/2,f(i)+kf2*dt/2,p(i)+kp2*dt/2,dth(i)+dkth2*dt/2,df(i)+dkf2*dt/2,dp(i)+dkp2*dt/2,I1,I3,hcm,M);
    %
     kth4=kth1+dkth3*dt;
    kf4=kf1+dkf3*dt;
    kp4=kp1+dkp3*dt;
    [dkth4,dkf4,dkp4]=peonza(th(i)+kth3*dt,f(i)+kf3*dt,p(i)+kp3*dt,dth(i)+dkth3*dt,df(i)+dkf3*dt,dp(i)+dkp3*dt,I1,I3,hcm,M);
    %
    th(i+1)=th(i)+(dt/6)*(kth1+2*kth2+2*kth3+kth4);
    f(i+1)=f(i)+(dt/6)*(kf1+2*kf2+2*kf3+kf4);
    p(i+1)=p(i)+(dt/6)*(kp1+2*kp2+2*kp3+kp4);
    %
        dth(i+1)=dth(i)+(dt/6)*(dkth1+2*dkth2+2*dkth3+dkth4);
    df(i+1)=df(i)+(dt/6)*(dkf1+2*dkf2+2*dkf3+dkf4);
    dp(i+1)=dp(i)+(dt/6)*(dkp1+2*dkp2+2*dkp3+dkp4);
    %
    [ddth(i+1),ddf(i+1),ddp(i+1)]=peonza(th(i+1),f(i+1),p(i+1),dth(i+1),df(i+1),dp(i+1),I1,I3,hcm,M);
end
te=linspace(0,2*pi,30);
r=linspace(0,Ra,30);
[T,R]=meshgrid(te,r);
XXc=R.*cos(T);
YYc=R.*sin(T);
ZZc=h*ones(length(r));
XX=R.*cos(T);
YY=R.*sin(T);
ZZ=h*R/Ra;
xa=0; ya=0; za=h*1.4; i=1; xb=Ra; yb=0; zb=h;
B=[cos(f(i)),-sin(f(i)),0;...
    sin(f(i)),cos(f(i)),0;...
    0,0,1];
D=[cos(p(i)),-sin(p(i)),0;...
    sin(p(i)),cos(p(i)),0;...
    0,0,1];
C=[1,0,0;...
    0,cos(th(i)),-sin(th(i));...
    0,sin(th(i)),cos(th(i))];
A=B*C*D;
V=A*[xa;ya;za];
xi=V(1); yi=V(2); zi=V(3);
%
V=A*[xb;yb;zb];
xib=V(1); yib=V(2); zib=V(3);
%%
n=1;
fig=figure(1);
pause(1.3)
i=N;
for i=1:N;
B=[cos(f(i)),-sin(f(i)),0;...
    sin(f(i)),cos(f(i)),0;...
    0,0,1];
D=[cos(p(i)),-sin(p(i)),0;...
    sin(p(i)),cos(p(i)),0;...
    0,0,1];
C=[1,0,0;...
    0,cos(th(i)),-sin(th(i));...
    0,sin(th(i)),cos(th(i))];
A=B*C*D;
Cero=A*[0;0;h];
Cerob=A*[0;0;h/10];
for j=1:30;
    for k=1:30;
        V=A*[XX(j,k);YY(j,k);ZZ(j,k)];
        V2=A*[XXc(j,k);YYc(j,k);ZZc(j,k)];
        XXv(j,k)=V(1);YYv(j,k)=V(2); ZZv(j,k)=V(3);
         XXd(j,k)=V2(1);YYd(j,k)=V2(2); ZZd(j,k)=V2(3);
    end
end
Va=A*[xa;ya;za];
xav=Va(1);yav=Va(2); zav=Va(3);
Va=A*[xb;yb;zb];
xbv=Va(1);ybv=Va(2); zbv=Va(3);
h1=surf(XXv,YYv,ZZv)
hold on
h2=surf(XXd,YYd,ZZd)
hold on
plot3([Cerob(1),xav],[Cerob(2),yav],[Cerob(3),zav],'-g','LineWidth',2)
hold on
plot3(xav,yav,zav,'.g','MarkerSize',5)
hold on
xi=[xi,xav];yi=[yi,yav];zi=[zi,zav];
plot3(xi(1:i),yi(1:i),zi(1:i),'-r','LineWidth',1.2)
hold on
%
plot3([Cero(1),xbv],[Cero(2),ybv],[Cero(3),zbv],'-g','LineWidth',2)
hold on
plot3(xbv,ybv,zbv,'.g','MarkerSize',5)
hold on
xib=[xib,xbv];yib=[yib,ybv];zib=[zib,zbv];
%plot3(xib(1:i),yib(1:i),zib(1:i),'-g')
camlight
axis([-2.7*Ra 2.7*Ra -2.7*Ra 2.7*Ra 0 1.7*h])
set(h1,'FaceColor','blue','FaceAlpha',1,'EdgeColor','none')
set(h2,'FaceColor','blue','FaceAlpha',1,'EdgeColor','none')
view(30,15)
hold off
colormap hsv
%set(gcf,'InvertHardCopy','off','Color','white');
daspect([1 1 1])
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
%axis off
frame=getframe(fig);
im=frame2im(frame); [imind,cm]=rgb2ind(im,256); imwrite(imind,cm,strcat(num2str(10000+n),'.png'));
n=n+1;
end
