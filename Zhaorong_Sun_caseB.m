clear; clc; close all;

filename = 'caseB_grid_battery_market_hourly.csv';
opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve';
data = readtable(filename, opts);

if ~isdatetime(data.timestamp)
    data.timestamp = datetime(data.timestamp, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
end

price_MWh = data.('day_ahead_price_gbp_per_mwh');
price_kWh = price_MWh / 1000;

T = height(data);
dt = 1;

Emax = 2000;
PchMax = 1000;
PdisMax = 1000;
eta_ch = 0.938;
eta_dis = 0.938;
SOC0 = 1000;
c_deg = 0.01;

nVar = 4 * T;
iCh  = 1:T;
iDis = T + (1:T);
iSOC = 2*T + (1:T);
iU   = 3*T + (1:T);
intcon = iU;

lb = zeros(nVar,1);
ub = inf(nVar,1);
ub(iCh) = PchMax;
ub(iDis) = PdisMax;
ub(iSOC) = Emax;
ub(iU) = 1;

Aeq = sparse(T + 1, nVar);
beq = zeros(T + 1, 1);

Aeq(1, iSOC(1)) = 1;
Aeq(1, iCh(1))  = -eta_ch * dt;
Aeq(1, iDis(1)) = dt / eta_dis;
beq(1) = SOC0;

for t = 2:T
    Aeq(t, iSOC(t))   = 1;
    Aeq(t, iSOC(t-1)) = -1;
    Aeq(t, iCh(t))    = -eta_ch * dt;
    Aeq(t, iDis(t))   = dt / eta_dis;
end

Aeq(T+1, iSOC(T)) = 1;
beq(T+1) = SOC0;

A = sparse(2*T, nVar);
b = zeros(2*T, 1);
for t = 1:T
    A(t,   iCh(t)) = 1;
    A(t,   iU(t))  = -PchMax;
    b(t) = 0;

    A(T+t, iDis(t)) = 1;
    A(T+t, iU(t))   = PdisMax;
    b(T+t) = PdisMax;
end

f_base = zeros(nVar, 1);
f_base(iCh)  = price_kWh * dt;
f_base(iDis) = -price_kWh * dt;

options = optimoptions('intlinprog', 'Display', 'final');
[x_base, fval_base, exitflag_base] = intlinprog(f_base, intcon, A, b, Aeq, beq, lb, ub, options);

if exitflag_base <= 0
    error('Base-case optimisation did not converge.');
end

Pch_base  = x_base(iCh);
Pdis_base = x_base(iDis);
SOC_base  = x_base(iSOC);
u_base    = x_base(iU);
netDispatch_base = Pdis_base - Pch_base;

profit_base = -fval_base;
throughput_base = sum((Pch_base + Pdis_base) * dt);

socBoundViolation_base = any(SOC_base < -1e-6 | SOC_base > Emax + 1e-6);
chargePowerViolation_base = any(Pch_base < -1e-6 | Pch_base > PchMax + 1e-6);
dischargePowerViolation_base = any(Pdis_base < -1e-6 | Pdis_base > PdisMax + 1e-6);
simultaneousViolation_base = any((Pch_base > 1e-6) & (Pdis_base > 1e-6));
terminalSOCViolation_base = abs(SOC_base(end) - SOC0) > 1e-6;

socCalc_base = zeros(T,1);
socCalc_base(1) = SOC0 + eta_ch * Pch_base(1) * dt - (Pdis_base(1) * dt) / eta_dis;
for t = 2:T
    socCalc_base(t) = SOC_base(t-1) + eta_ch * Pch_base(t) * dt - (Pdis_base(t) * dt) / eta_dis;
end
socEquationMismatch_base = max(abs(socCalc_base - SOC_base)) > 1e-5;

summaryBase = table(...
    ["Total profit (GBP)"; ...
     "Total throughput (kWh)"; ...
     "Initial SOC (kWh)"; ...
     "Final SOC (kWh)"; ...
     "Minimum SOC (kWh)"; ...
     "Maximum SOC (kWh)"; ...
     "SOC bound violation"; ...
     "Charge power violation"; ...
     "Discharge power violation"; ...
     "Simultaneous charge/discharge"; ...
     "Terminal SOC violation"; ...
     "SOC equation mismatch"], ...
    [profit_base; throughput_base; SOC0; SOC_base(end); min(SOC_base); max(SOC_base); ...
     double(socBoundViolation_base); double(chargePowerViolation_base); ...
     double(dischargePowerViolation_base); double(simultaneousViolation_base); ...
     double(terminalSOCViolation_base); double(socEquationMismatch_base)], ...
    'VariableNames', {'Metric','Value'});

writetable(summaryBase, 'caseB_summary_results.csv');

dispatchBase = table(data.timestamp, price_MWh, price_kWh, Pch_base, Pdis_base, netDispatch_base, SOC_base, u_base, ...
    'VariableNames', {'timestamp','price_MWh','price_kWh','Pch_kW','Pdis_kW','netDispatch_kW','SOC_kWh','u'});
writetable(dispatchBase, 'caseB_dispatch_results.csv');

f_ext = zeros(nVar, 1);
f_ext(iCh)  = (price_kWh + c_deg) * dt;
f_ext(iDis) = (c_deg - price_kWh) * dt;

[x_ext, fval_ext, exitflag_ext] = intlinprog(f_ext, intcon, A, b, Aeq, beq, lb, ub, options);

if exitflag_ext <= 0
    error('Degradation-case optimisation did not converge.');
end

Pch_ext  = x_ext(iCh);
Pdis_ext = x_ext(iDis);
SOC_ext  = x_ext(iSOC);
u_ext    = x_ext(iU);
netDispatch_ext = Pdis_ext - Pch_ext;

throughput_ext = sum((Pch_ext + Pdis_ext) * dt);
grossProfit_ext = sum(price_kWh .* (Pdis_ext - Pch_ext) * dt);
degradationCost_ext = c_deg * throughput_ext;
netProfit_ext = grossProfit_ext - degradationCost_ext;

socBoundViolation_ext = any(SOC_ext < -1e-6 | SOC_ext > Emax + 1e-6);
chargePowerViolation_ext = any(Pch_ext < -1e-6 | Pch_ext > PchMax + 1e-6);
dischargePowerViolation_ext = any(Pdis_ext < -1e-6 | Pdis_ext > PdisMax + 1e-6);
simultaneousViolation_ext = any((Pch_ext > 1e-6) & (Pdis_ext > 1e-6));
terminalSOCViolation_ext = abs(SOC_ext(end) - SOC0) > 1e-6;

socCalc_ext = zeros(T,1);
socCalc_ext(1) = SOC0 + eta_ch * Pch_ext(1) * dt - (Pdis_ext(1) * dt) / eta_dis;
for t = 2:T
    socCalc_ext(t) = SOC_ext(t-1) + eta_ch * Pch_ext(t) * dt - (Pdis_ext(t) * dt) / eta_dis;
end
socEquationMismatch_ext = max(abs(socCalc_ext - SOC_ext)) > 1e-5;

summaryExt = table(...
    ["Degradation penalty (GBP per kWh)"; ...
     "Gross trading profit (GBP)"; ...
     "Degradation cost (GBP)"; ...
     "Net profit (GBP)"; ...
     "Total throughput (kWh)"; ...
     "Initial SOC (kWh)"; ...
     "Final SOC (kWh)"; ...
     "Minimum SOC (kWh)"; ...
     "Maximum SOC (kWh)"; ...
     "SOC bound violation"; ...
     "Charge power violation"; ...
     "Discharge power violation"; ...
     "Simultaneous charge/discharge"; ...
     "Terminal SOC violation"; ...
     "SOC equation mismatch"], ...
    [c_deg; grossProfit_ext; degradationCost_ext; netProfit_ext; throughput_ext; SOC0; SOC_ext(end); ...
     min(SOC_ext); max(SOC_ext); double(socBoundViolation_ext); double(chargePowerViolation_ext); ...
     double(dischargePowerViolation_ext); double(simultaneousViolation_ext); ...
     double(terminalSOCViolation_ext); double(socEquationMismatch_ext)], ...
    'VariableNames', {'Metric','Value'});

writetable(summaryExt, 'caseB_extension_summary_results.csv');

dispatchExt = table(data.timestamp, price_MWh, price_kWh, Pch_ext, Pdis_ext, netDispatch_ext, SOC_ext, u_ext, ...
    'VariableNames', {'timestamp','price_MWh','price_kWh','Pch_kW','Pdis_kW','netDispatch_kW','SOC_kWh','u'});
writetable(dispatchExt, 'caseB_extension_dispatch_results.csv');

comparisonTable = table(...
    ["Total profit (GBP)"; "Total throughput (kWh)"; "Final SOC (kWh)"], ...
    [profit_base; throughput_base; SOC_base(end)], ...
    [netProfit_ext; throughput_ext; SOC_ext(end)], ...
    'VariableNames', {'Metric','Base_case','Degradation_case'});
writetable(comparisonTable, 'caseB_base_vs_extension_comparison.csv');

figure(1);
yyaxis left
plot(data.timestamp, price_MWh, 'LineWidth', 1);
ylabel('Day-ahead price (GBP/MWh)');
yyaxis right
stairs(data.timestamp, netDispatch_base, 'LineWidth', 1);
ylabel('Net dispatch (kW)');
xlabel('Time');
title('Day-ahead price and battery dispatch');
grid on;
saveas(gcf, 'figure1_price_dispatch.png');

figure(2);
plot(data.timestamp, SOC_base, 'LineWidth', 1);
xlabel('Time');
ylabel('SOC (kWh)');
title('Battery SOC profile');
grid on;
saveas(gcf, 'figure2_soc.png');

idx = data.timestamp >= datetime(2025,6,1) & data.timestamp < datetime(2025,6,8);

figure(3);
yyaxis left
plot(data.timestamp(idx), price_MWh(idx), 'LineWidth', 1);
ylabel('Day-ahead price (GBP/MWh)');
yyaxis right
stairs(data.timestamp(idx), netDispatch_base(idx), 'LineWidth', 1);
ylabel('Net dispatch (kW)');
xlabel('Time');
title('Day-ahead price and battery dispatch (7-day window)');
grid on;
saveas(gcf, 'figure1_price_dispatch_7days.png');

figure(4);
plot(data.timestamp(idx), SOC_base(idx), 'LineWidth', 1.2);
xlabel('Time');
ylabel('SOC (kWh)');
title('Battery SOC profile (7-day window)');
grid on;
saveas(gcf, 'figure2_soc_7days.png');

figure(5);
yyaxis left
plot(data.timestamp(idx), price_MWh(idx), 'LineWidth', 1);
ylabel('Day-ahead price (GBP/MWh)');
yyaxis right
stairs(data.timestamp(idx), netDispatch_ext(idx), 'LineWidth', 1);
ylabel('Net dispatch (kW)');
xlabel('Time');
title('Day-ahead price and battery dispatch with degradation (7-day window)');
grid on;
saveas(gcf, 'figure3_price_dispatch_7days_degradation.png');

figure(6);
plot(data.timestamp(idx), SOC_ext(idx), 'LineWidth', 1.2);
xlabel('Time');
ylabel('SOC (kWh)');
title('Battery SOC profile with degradation (7-day window)');
grid on;
saveas(gcf, 'figure4_soc_7days_degradation.png');

save('caseB_clean_full_workspace.mat');
