%--------------------------------------------------------------------------
%Used to fix some parameters and let the others vary (INDMAP) before
%solving ODE
%--------------------------------------------------------------------------

function [rout,J,CV] = DriverBasic_plot(pars,data,ptnb)

[rout,J,CV] = CVmodel(pars,data);

N   = length(data.td_per);
tdc = data.td_per;
pla = CV.plaS(end-N+1:end);
plv = CV.plvS(end-N+1:end);
psa = CV.psaS(end-N+1:end);
psv = CV.psvS(end-N+1:end);
pra = CV.praS(end-N+1:end);
prad= CV.pRAdS(end-N+1:end);
prv = CV.prvS(end-N+1:end);
prvd= CV.pRVdS(end-N+1:end);
ppa = CV.ppaS(end-N+1:end);
ppad= CV.pPAdS(end-N+1:end);
ppv = CV.ppvS(end-N+1:end);
CO  = CV.COS(end);

qava = CV.qava(end-N+1:end);
qmva = CV.qmva(end-N+1:end);
qtva = CV.qtva(end-N+1:end);
qsv  = CV.qsv(end-N+1:end);
qpva = CV.qpva(end-N+1:end);

Vrv = CV.VrvS(end-N+1:end);
Vra = CV.VraS(end-N+1:end);
Vlv = CV.VlvS(end-N+1:end);
Vla = CV.VlaS(end-N+1:end);


states.p = [pla; plv; psa; psv; pra; prv; ppa; ppv];
states.q = [qpva; qmva; qtva; qsv; qpva];
states.V = [Vla; Vlv; Vra; Vrv];


figure(ptnb+200);hold on;
%%% pLA
subplot(3,3,1); hold on
h = plot(tdc,pla,'b');
set(h,'LineWidth',2);
set(gca,'FontSize',20);
ylabel('pla (mmHg)');

%%% pLV
subplot(3,3,2);hold on
h = plot(tdc,plv,'b');
set(h,'LineWidth',2);
set(gca,'FontSize',20);
ylabel('plv (mmHg)');
%     title(st);

%%% pSA
maxpSA = data.pSAM*ones(size(tdc));
minpSA = data.pSAm*ones(size(tdc));
subplot(3,3,3); hold on
h = plot(tdc,psa,'b',tdc,maxpSA,'--k',tdc,minpSA,'--k');
set(h,'Linewidth',2);
set(gca,'FontSize',20);
ylabel('psa (mmHg)');

%%% pSV
subplot(3,3,4); hold on
h = plot(tdc,psv,'b');
set(h,'LineWidth',2);
set(gca,'FontSize',20);
ylabel('psv (mmHg)');

%%% pRA
Nra = round(2*length(data.pRA)/3);
maxpRA = max(data.pRA(Nra:end))*ones(size(tdc));
minpRA = min(data.pRA(1:round(length(data.pRA)/3)));
subplot(3,3,5); hold on
h = plot(tdc,pra,'b',tdc,prad,'k',tdc,maxpRA,'--k',tdc,minpRA,'--k');
set(h,'linewidth',2);
set(gca,'FontSize',20);
ylabel('pra (mmHg)');

%%% pRV
maxpRV=max(data.pRV)*ones(size(tdc));
minpRV=min(data.pRV)*ones(size(tdc));
subplot(3,3,6); hold on
h = plot(tdc,prv,'b',tdc,prvd,'k',tdc, maxpRV,'--k',tdc,minpRV, '--k');
set(h,'LineWidth',2);
set(gca,'FontSize',20);
ylabel('prv (mmHg)');

%%% pPA
maxpPA = max(data.pPA)*ones(size(tdc));
minpPA = min(data.pPA)*ones(size(tdc));
subplot(3,3,7);hold on;
h = plot(tdc,ppa,'b',tdc,circshift(ppad,-2),'k',tdc,maxpPA,'--k',tdc,minpPA,'--k');
set(h,'LineWidth',2);
set(gca,'FontSize',20);
ylabel('ppa (mmHg)');
xlabel('time (s)');

%%% pPV
minppv = min(data.pPW)*ones(size(tdc));
subplot(3,3,8);hold on
h = plot(tdc,ppv,'b',tdc,minppv,'--k');
set(h,'Linewidth',2);
set(gca,'FontSize',20);
ylabel('ppv (mmHg)');
xlabel('time (s)');

%%% CO
COt = CO*ones(size(tdc));
COdt = data.CO.*ones(size(tdc));
subplot(3,3,9); hold on
h = plot(tdc,COt,'b',tdc,COdt,'--k');
set(h,'linewidth',2);
set(gca,'FontSize',20);
ylabel('CO (mL/sec)');
xlabel('time (s)');

figure(ptnb+20);hold on;
subplot(1,2,1);hold on;
h=plot(Vrv,prv,'b',Vlv,plv,'r');
set(h,'linewidth',3);
set(gca,'fontsize',20);
grid on;
xlabel('Ventricle Volume (mL)');
ylabel('Ventricle Pressure (mmHg)');
legend('RV','LV');
subplot(1,2,2);hold on;
h=plot(Vra,pra,'b',Vla,pla,'r');
set(h,'linewidth',3);
set(gca,'fontsize',20);
grid on;
legend('RA','LA','Location','northwest');
xlabel('Atrium Volume (mL)');
ylabel('Atrium Pressure (mmHg)');
grid on;

end
