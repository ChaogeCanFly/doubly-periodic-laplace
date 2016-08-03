function F=Video(X,Y,H,Frames,Xps)
% This program displays a video of the results.  X and Y are the
% corresponding results from Membrane.m.  The parallel walls are y=H and
% y=-H.  Frames is the number of frames per second.
X=X(~cellfun('isempty',X));
disp([num2str(length(X)),' frames'])
Y=Y(~cellfun('isempty',Y));
n=length(X);
M=length(X{1});
Xc=cell(1,M);
Yc=cell(1,M);
t=linspace(0,2*pi,129*20)';
for k=1:M
    for j=1:n
        Xc{k}=[Xc{k},X{j}{k}];
        Yc{k}=[Yc{k},Y{j}{k}];
    end
end
XF=cell(1,M);
YF=cell(1,M);

P=linspace(0,2*pi,51)';
for k=1:M
    XF{k}=It(Xc{k},P);
    YF{k}=It(Yc{k},P);
end

k=0;
m=1;
[U,D]=MakeUD(2000);
Uout=U.Z(t);
Dout=D.Z(t);
x1=real(Uout);
y1=imag(Uout);
x2=real(Dout);
y2=imag(Dout);
figure('units','normalized','outerposition',[0 0 1 1])
if nargout==0
    %h=figure(1);
    %set(h, 'Position', [150 400 1000/2 round(1000*H/pi)])
    %     set(gcf,'PaperPositionMode','auto')
    while k<n
        % while k<n
        %k=k+1;
        k=mod(k,n-1)+1;
        if nargin>4
            if k<2
                k2=2;
            else
                k2=k;
            end
            for j=1:length(Xps{1})
            %quiver(Xps{j}(:,1),Xps{j}(:,2),Vps{k2-1}{j}(:,1)/4,Vps{k2-1}{j}(:,2)/4,0,'color',[0 0.5 0],'ShowArrowHead','off')
            scatter(Xps{k2-1}{j}(:,1),Xps{k2-1}{j}(:,2))
            hold on
            end
            plot(x1,y1)
        else
            plot(x1,y1)
            hold on
        end
        plot(x2,y2)
        for j=1:M
            plot(XF{j}(:,k),YF{j}(:,k),'r-')
            plot(XF{j}(:,k)-2*pi,YF{j}(:,k),'r-')
            plot(XF{j}(:,k)+2*pi,YF{j}(:,k),'r-')
            %plot(XF{j}(:,k)+4*pi,YF{j}(:,k),'r-')
        end
        %title([num2str(round(100*(k-1)/(n-1))),'%'])
        title(num2str(k))
        axis equal
        axis([0, 2*pi -H H]) % We assume the membrane is in this box.
        hold off
        axis off
        %         if k<10
        %             print(h,'-dpng','-r500',['image000',num2str(k),'.png'])
        %         elseif k<100
        %             print(h,'-dpng','-r500',['image00',num2str(k),'.png'])
        %         elseif k<1000
        %             print(h,'-dpng','-r500',['image0',num2str(k),'.png'])
        %         else
        %             print(h,'-dpng','-r500',['image',num2str(k),'.png'])
        %         end
        pause(1/Frames)
        if k==1
            pause(1)
        end
    end
else
    while k<n-1
        k=k+1;
        if nargin>4
            if k==1
                k2=2;
            else
                k2=k;
            end
            for j=1:length(Xps{1})
            %quiver(Xps{j}(:,1),Xps{j}(:,2),Vps{k2-1}{j}(:,1)/4,Vps{k2-1}{j}(:,2)/4,0,'color',[0 0.5 0],'ShowArrowHead','off')
            scatter(Xps{k2-1}{j}(:,1),Xps{k2-1}{j}(:,2))
            hold on
            end
            plot(x1,y1)
        else
            plot(x1,y1)
            hold on
        end
        plot(x2,y2)
        for j=1:M
            plot(XF{j}(:,k),YF{j}(:,k),'r-')
            plot(XF{j}(:,k)-2*pi,YF{j}(:,k),'r-')
            plot(XF{j}(:,k)+2*pi,YF{j}(:,k),'r-')
            %plot(XF{j}(:,k)+4*pi,YF{j}(:,k),'r-')
        end
        %title([num2str(round(100*(k-1)/(n-1))),'%'])
        title(num2str(k))
        axis equal
        axis([0, 2*pi -H H]) % We assume the membrane is in this box.
        hold off
        axis off
        pause(1/Frames)
        h=figure(1);
        F(m)=getframe(h);
        m=m+1;
        if k==1
            pause(1)
        end
    end
end
end