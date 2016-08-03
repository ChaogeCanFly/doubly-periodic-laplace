function [z]=Lap_DLP_Close_Apply_demo(NPt,s,t,f)

M=length(s.len);

if isfield(NPt,'s')
    z.sub.s=zeros(length(s.x),1);
    z.add.s=z.sub.s;
    z.sub.sx=z.sub.s;
    z.add.sx=z.sub.s;
    z.sub.sy=z.sub.s;
    z.add.sy=z.sub.s;
    ls=0;
    for k=1:M
        if ~isempty(NPt.s{k})
            ss.x=s.x(ls+1:ls+s.len(k));
            ss.t=s.t(ls+1:ls+s.len(k));
            ss.xp=s.xp(ls+1:ls+s.len(k));
            ss.sp=s.sp(ls+1:ls+s.len(k));
            ss.w=s.w(ls+1:ls+s.len(k));
            ss.nx=s.nx(ls+1:ls+s.len(k));
            ss.cw=s.cw(ls+1:ls+s.len(k));
            xx.x=s.x(NPt.s{k});
            ff=f(ls+1:ls+s.len(k));
            %in=inpoly([real(xx),imag(xx)],[real(ss.x),imag(ss.x)]);
            in=inpolygon(real(xx.x),imag(xx.x),real(ss.x),imag(ss.x));
            zz=zeros(length(xx.x),1);
            zzx=zz;
            zzy=zz;
            sin=sum(in);
            if sin>0
                disp('i1')
                [u,ux,uy]=lapDevalclose(xx.x(in),ss,ff,'i');
                zz(in)=u;
                zzx(in)=ux;
                zzy(in)=uy;
                %zzn(in)=real(xx.nx(in)).*ux+imag(xx.nx(in)).*uy;
            end
            if sin<length(in)
                disp('e1')
                [u,ux,uy]=lapDevalclose(xx.x(~in),ss,ff,'e');
                zz(~in)=u;
                zzx(~in)=ux;
                zzy(~in)=uy;
                %zzn(~in)=real(xx.nx(~in)).*ux+imag(xx.nx(~in)).*uy;
            end
            DLP=LapDLPmatrix(xx,ss,0)*ff;
            xx.nx=0*xx.x+1;
            [~,DLPx]=LapDLPmatrix(xx,ss,0);
            DLPx=DLPx*ff;
            xx.nx=0*xx.x+1i;
            [~,DLPy]=LapDLPmatrix(xx,ss,0);
            DLPy=DLPy*ff;
            z.sub.s(NPt.s{k})=z.sub.s(NPt.s{k})+DLP;
            z.add.s(NPt.s{k})=z.add.s(NPt.s{k})+zz;
            z.sub.sx(NPt.s{k})=z.sub.sx(NPt.s{k})+DLPx;
            z.add.sx(NPt.s{k})=z.add.sx(NPt.s{k})+zzx;
            z.sub.sy(NPt.s{k})=z.sub.sy(NPt.s{k})+DLPy;
            z.add.sy(NPt.s{k})=z.add.sy(NPt.s{k})+zzy;
        end
        ls=ls+s.len(k);
    end
end

if isfield(NPt,'t')
    z.sub.t=zeros(length(t.x),1);
    z.add.t=z.sub.t;
    z.sub.tx=z.sub.t;
    z.add.tx=z.sub.t;
    z.sub.ty=z.sub.t;
    z.add.ty=z.sub.t;
    ls=0;
    for k=1:M
        if ~isempty(NPt.t{k})
            ss.x=s.x(ls+1:ls+s.len(k));
            ss.t=s.t(ls+1:ls+s.len(k));
            ss.xp=s.xp(ls+1:ls+s.len(k));
            ss.sp=s.sp(ls+1:ls+s.len(k));
            ss.w=s.w(ls+1:ls+s.len(k));
            ss.nx=s.nx(ls+1:ls+s.len(k));
            ss.cw=s.cw(ls+1:ls+s.len(k));
            xx.x=t.x(NPt.t{k});
            ff=f(ls+1:ls+s.len(k));
            %in=inpoly([real(xx),imag(xx)],[real(ss.x),imag(ss.x)]);
            in=inpolygon(real(xx.x),imag(xx.x),real(ss.x),imag(ss.x));
            zz=zeros(length(xx.x),1);
            zzx=zz;
            zzy=zz;
            sin=sum(in);
            if sin>0
                disp('i')
                [u,ux,uy]=lapDevalclose(xx.x(in),ss,ff,'i');
                zz(in)=u;
                zzx(in)=ux;
                zzy(in)=uy;
            end
            if sin<length(in)
                disp('e')
                [u,ux,uy]=lapDevalclose(xx.x(~in),ss,ff,'e');
                zz(~in)=u;
                zzx(~in)=ux;
                zzy(~in)=uy;
            end
            DLP=LapDLPmatrix(xx,ss,0)*ff;
            xx.nx=0*xx.x+1;
            [~,DLPx]=LapDLPmatrix(xx,ss,0);
            DLPx=DLPx*ff;
            xx.nx=0*xx.x+1i;
            [~,DLPy]=LapDLPmatrix(xx,ss,0);
            DLPy=DLPy*ff;
            z.sub.t(NPt.t{k})=z.sub.t(NPt.t{k})+DLP;
            z.add.t(NPt.t{k})=z.add.t(NPt.t{k})+zz;
            z.sub.tx(NPt.t{k})=z.sub.tx(NPt.t{k})+DLPx;
            z.add.tx(NPt.t{k})=z.add.tx(NPt.t{k})+zzx;
            z.sub.ty(NPt.t{k})=z.sub.ty(NPt.t{k})+DLPy;
            z.add.ty(NPt.t{k})=z.add.ty(NPt.t{k})+zzy;
        end
        ls=ls+s.len(k);
    end
end