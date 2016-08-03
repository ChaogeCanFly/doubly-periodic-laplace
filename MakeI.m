function [I,N]=MakeI(clup,M,n)

%[I,N]=GridI;

[I,N]=FastRandI(clup,M,n);
% n=128;
% t=linspace(0,2*pi,n+1);
% t=t(2:end);
% 
% % Islands are defined here
% I{1}.x=cos(t)-1+1i*(sin(t)/2+cos(2*t)/4-1/2);
% N(1)=n;
% I{2}.x=cos(t)/2-1+1i*sin(t)/4+0.2*1i+1.4*1i+1;
% N(2)=n;
% I{3}.x=cos(t)/4+1i*sin(t)+1-1i;
% N(3)=n;
% I{4}.x=cos(t)/4-1+1i*sin(t)/4+0.2*1i+0.4*1i-0.5;
% N(4)=n;
% I{5}.x=cos(t)/6-1+1i*sin(t)/8+0.1*1i;
% N(5)=n;
% I{6}.x=cos(t)/2+1i*sin(t)/4-2-2.8*1i;
% N(6)=n;
% I{7}.x=cos(t)/4+1i*sin(t)/3+1.5-2.8*1i;
% N(7)=n;
% I{8}.x=cos(t)/8+1i*sin(t)/7+1+2.8*1i;
% N(8)=n;
% I{9}.x=cos(t)/8+1i*sin(t)/7+0.5+0.5*1i;
% N(9)=n;
% I{10}.x=cos(t)/8+1i*sin(t)/5-2+2*1i;
% N(10)=n;
% I{11}.x=cos(t)/8+1i*sin(t)/8-1-2*1i;
% N(11)=n;
% I{12}.x=cos(t)/8+1i*sin(t)/8+2.5*1i;
% N(12)=n;
% I{13}.x=cos(t)/4+1i*sin(t)/4-2.5*1i;
% N(13)=n;
% I{14}.x=cos(t)/4+1i*sin(t)+2+1.5i;
% N(14)=n;
% I{15}.x=cos(t)/8+1i*sin(t)/8+0.2+3*1i;
% N(15)=n;
% I{16}.x=cos(t)/8+1i*sin(t)/8+2-2*1i;
% N(16)=n;
% I{17}.x=cos(t)/8+1i*sin(t)/8+0.9*1i;
% N(17)=n;
% I{18}.x=cos(t)/8+1i*sin(t)/8+1+0.2*1i;
% N(18)=n;
% I{19}.x=cos(t)/8+1i*sin(t)/8-1.5*1i;
% N(19)=n;
% I{20}.x=cos(t)/8+1i*sin(t)/8+1-2.3*1i;
% N(20)=n;
% I{21}.x=cos(t)/8+1i*sin(t)/8+1.8+2.9*1i;
% N(21)=n;
% I{22}.x=cos(t)/4+1i*sin(t)/4+1.75-0.2*1i;
% N(22)=n;
% I{23}.x=cos(t)/8+1i*sin(t)/8+2.2+0.2*1i;
% N(23)=n;
% I{24}.x=cos(t)/8+1i*sin(t)/8+0.8-2.8*1i;
% N(24)=n;
% I{25}.x=cos(t)/8+1i*sin(t)/8+1.4+2.4*1i;
% N(25)=n;
% I{26}.x=cos(t)/8+1i*sin(t)/8-2-1.5*1i;
% N(26)=n;
% I{27}.x=cos(t)/8+1i*sin(t)/8-1.8+1.2*1i;
% N(27)=n;
% I{28}.x=cos(t)/4+1i*sin(t)/4-1.2+2.5*1i;
% N(28)=n;
% I{29}.x=cos(t)/2+1i*sin(t)/2+2.5-1*1i;
% N(29)=n;
% I{30}.x=cos(t)/4+1i*sin(t)/4-2-2.2*1i;
% N(30)=n;
% I{31}.x=cos(t)/4+1i*sin(t)/4-1.5+1.8*1i;
% N(31)=n;
% I{32}.x=cos(t)/4+1i*sin(t)/4+2.5-2.7*1i;
% N(32)=n;
% I{33}.x=cos(t)/8+1i*sin(t)/8-1.2-1.5*1i;
% N(33)=n;
% I{34}.x=cos(t)/8+1i*sin(t)/8-0.2+0.5*1i;
% N(34)=n;
% I{35}.x=cos(t)/8+1i*sin(t)/8-0.8-2.7*1i;
% N(35)=n;
% I{36}.x=cos(t)/8+1i*sin(t)/8-0.5-1.9*1i;
% N(36)=n;
% I{37}.x=cos(t)/8+1i*sin(t)/8+0.2*1i;
% N(37)=n;
% I{38}.x=cos(t)/4+1i*sin(t)/4+2.7+2.7*1i;
% N(38)=n;
% I{39}.x=cos(t)/4+1i*sin(t)/4-1.1+1.2*1i;
% N(39)=n;
% I{40}.x=cos(t)/4+1i*sin(t)/4-0.6-2.3*1i;
% N(40)=n;
% I{41}.x=cos(t)/8+1i*sin(t)/8-0.6+2.8*1i;
% N(41)=n;
% I{42}.x=cos(t)/8+1i*sin(t)/8-1.6+2.9*1i;
% N(42)=n;
% I{43}.x=cos(t)/8+1i*sin(t)/8-0.7+2.1*1i;
% N(43)=n;
% I{44}.x=cos(t)/8+1i*sin(t)/8-0.6-1.5*1i;
% N(44)=n;
% I{45}.x=cos(t)/2+1i*sin(t)/2+1.1+1.1*1i;
% N(45)=n;
% I{46}.x=cos(t)/4+1i*sin(t)/4+2.7-1.9*1i;
% N(46)=n;
% I{47}.x=cos(t)/4+1i*sin(t)/4+0.8+2.2*1i;
% N(47)=n;
% I{48}.x=cos(t)/4+1i*sin(t)/4+2.8+1*1i;
% N(48)=n;
% I{49}.x=cos(t)/4+1i*sin(t)/4+0.3-0.8*1i;
% N(49)=n;
% I{50}.x=cos(t)/4+1i*sin(t)/4+0.5-1.8*1i;
% N(50)=n;

% % %Plots islands for convenience
% figure(2)
% for k=1:length(N)
%    for i=-1:1
%       for j=-1:1
%           fill(real(I{k}.x)+2*pi*i,imag(I{k}.x)+2*pi*j,[0.7 0.7 0.7])
%           hold on
%       end
%    end
% end
% axis equal
% axis([-pi pi -pi pi])
% pause(0.1)