function [YVMDECDF] = VMDECDFbiasCorrection(Obsij,Satij,barrierPos,dayInterval)
    Ycal=Obsij(1:barrierPos);                       Yval=Obsij(1+barrierPos:end);
    Xcal=Satij(1:barrierPos);                       Xval=Satij(1+barrierPos:end);

    [Yimf,~,~] = emd(Ycal);
    nimf=1*size(Yimf,2);
    if size(Yimf,2) ~= 0
        [Ximf,Xresidual] = vmd(Xcal,'NumIMF',nimf,'MaxIterations',1000);
        [Yimf,Yresidual] = vmd(Ycal,'NumIMF',nimf,'MaxIterations',1000);
    else
        YVMDECDF= biasCorrection1D('ECDF',Obsij,Satij,barrierPos,[]);
        return
    end
    Ycalnew=cat(2,Yimf,Yresidual);      Xcalnew=cat(2,Ximf,Xresidual);
    clear Yimf Yresidual Ximf Xresidual
    Ximfval=nan(length(Yval),nimf);            Xresidualval=nan(length(Yval),1);
    Yimfval=nan(length(Yval),nimf);            Yresidualval=nan(length(Yval),1);
    
    if dayInterval == 0
        [Xtempimf,Xtempresidual] = vmd([Xcal;Xval],'NumIMF',nimf,'MaxIterations',1000);
        [Ytempimf,Ytempresidual] = vmd([Ycal;Yval],'NumIMF',nimf,'MaxIterations',1000);
        Ximfval=Xtempimf(barrierPos+1:end,:);
        Xresidualval=Xtempresidual(barrierPos+1:end,:);
        Yimfval=Ytempimf(barrierPos+1:end,:);
        Yresidualval=Ytempresidual(barrierPos+1:end,:);
        clear Xtempimf Xtempresidual Ytempimf Ytempresidual

    elseif dayInterval == 1
        for i=1:length(Yval)
            [Xtempimf,Xtempresidual] = vmd([Xcal;Xval(i)],'NumIMF',nimf,'MaxIterations',1000);
            [Ytempimf,Ytempresidual] = vmd([Ycal;Yval(i)],'NumIMF',nimf,'MaxIterations',1000);
            Ximfval(i,:)=Xtempimf(end,:);
            Xresidualval(i,:)=Xtempresidual(end,:);
            Yimfval(i,:)=Ytempimf(end,:);
            Yresidualval(i,:)=Ytempresidual(end,:);
            clear Xtempimf Xtempresidual Ytempimf Ytempresidual
        end
        
    else
        for i=1:dayInterval:length(Yval)
            epos=i+dayInterval-1;
            if epos > length(Yval)
                epos=length(Yval);
            end
            [Xtempimf,Xtempresidual] = vmd([Xcal;Xval(1:epos)],'NumIMF',nimf,'MaxIterations',1000);
            [Ytempimf,Ytempresidual] = vmd([Ycal;Yval(1:epos)],'NumIMF',nimf,'MaxIterations',1000);
            Ximfval(i:epos,:)=Xtempimf(end-(epos-i):end,:);
            Xresidualval(i:epos,:)=Xtempresidual(end-(epos-i):end,:);
            Yimfval(i:epos,:)=Ytempimf(end-(epos-i):end,:);
            Yresidualval(i:epos,:)=Ytempresidual(end-(epos-i):end,:);
            clear Xtempimf Xtempresidual Ytempimf Ytempresidual
        end
    end
    Yvalnew=cat(2,Yimfval,Yresidualval);      Xvalnew=cat(2,Ximfval,Xresidualval);
    clear Yimfval Yresidualval  Ximfval Xresidualval
    Ynew=cat(1,Ycalnew,Yvalnew);
    Xnew=cat(1,Xcalnew,Xvalnew);
    clear Ycalnew Yvalnew Xcalnew Xvalnew
    YECDFimf=Ynew*nan;
    for i=1:size(Ynew,2)
        [ YECDFimf(:,i)] = biasCorrection1D('ECDF',Ynew(:,i),Xnew(:,i),barrierPos,[]);
    end
    YVMDECDF=sum(YECDFimf,2);
    clear YECDFimf Ynew Xnew
    YVMDECDF(YVMDECDF<0)=0;

end