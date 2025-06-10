%% 1. LS ECDF LOCI QM RETree SVM
tic
clear BCName SatDataName
Selected_SatDataName={'IMERGEatIMD'}; %SatDataName={'IMERGEatIMD';'IMERGLatIMD';'IMERGFatIMD'};
DataName={'IMERGE'};
addpath(genpath('C:\Users\Shushobhit\Dropbox\schaudhary_codes'))
[~,~,n3]=size(IMD);
if strcmp(case_season,'JJAS')
   calValBarrier=round(n3*13/19);
elseif strcmp(case_season,'annual')
   TempDateCharMatrix = datestr(datenum(analysis_start_date):datenum(analysis_end_date),'yyyy-mm-dd');
   calValBarrier=1+diff([datenum(analysis_start_date) datenum('2013-12-31')]);
end
BCTech={'ECDF'};
BCName={'ECDF'};
BCVMDECDF=nan*IMD;
MySat=eval(char(Selected_SatDataName));
for r=1:size(IMD,1)
    for c=1:size(IMD,2)
        Y=(squeeze(IMD(r,c,:)));                                        X=(squeeze(MySat(r,c,:)));
%                 Y=flipud(Y);                                        X=flipud(X);
        if sum(~isnan(Y))>10
                  disp(['Row =' num2str(r) ' & Coloumn =' num2str(c)])

            tic
            BCVMDECDF(r,c,:) = VMDECDFbiasCorrection(Y,X,calValBarrier,0);
           toc
            % plotObsSatIMFs(Yorig,YVMDECDF,Ynew,YECDFimf);
            % figure; plot(YVMDECDF) ;hold on; plot(Y,'--r');
            % mse(YVMDECDF,Y)
        end
    end
end
save([char(DataName) '_' case_season '_BCVMDECDF'],'BCVMDECDF','calValBarrier','-v7.3')