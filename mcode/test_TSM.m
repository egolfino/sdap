function test_TSM
rawDataFN='Y:\CS_2004\SAS\AS18\end\rep4\trial-18-1.mat';

load(rawDataFN);    % gives data

MexIO('reset');
% data.params.rmsThresh=data.params.rmsThresh/2;
fmts=testTSM(data,data.params);

%%

figure;
subplot(2,1,1);
plot(data.fmts(:,1:2));
subplot(2,1,2);
% plot(data.rms(:,1));
% hold on;
% xs=get(gca,'XLim');
% plot(xs,repmat(data.params.rmsThresh,1,2),'k-');
plot(data.rms(:,2)./data.rms(:,3));
hold on;
xs=get(gca,'XLim');
plot(xs,repmat(data.params.rmsRatioThresh,1,2),'k-');


figure;
subplot(2,1,1);
plot(fmts);
return