function batch_fyn_gf()
gf_con = [80.65, 16.13, 8.065, 4.033, 1.515, 0.76, 0.1515]; % growth factor concentration in nM
% legend - growth factor concentration in ng/ml
leg = {'500 ng/ml', '100 ng/ml', '50 ng/ml', '25 ng/ml', '10 ng/ml', '5 ng/ml', '1 ng/ml'};
line_type = {'r-', 'g-', 'b-', 'k-', 'ro', 'go','bo', 'ko'};
n = length( gf_con);
figure; hold on;
for i = 1:n,
    [t, y] = fyn_gf_test('model', 'egfr_huang', 'gf_1', gf_con(i), 'show_figure', 0); % 50 ng/ml
    plot(t/60, y, line_type{i});
    clear t y;
    xlabel('Time (min)'); ylabel('Fyn Activity (nM)');
end;
legend('500 ng/ml', '100 ng/ml', '50 ng/ml', '25 ng/ml', '10 ng/ml', '5 ng/ml', '1 ng/ml');
return;