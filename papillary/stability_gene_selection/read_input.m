function rank_scr_cho_supl

fpqm =  csvread('../data/fpqm .csv', 1,1);
fpqm_log =  csvread('../data/fpqm_log .csv', 1,1);
nt = csvread('../data/nt .csv',1,1);
vs = csvread('../data/vs .csv',1,1);
labels = ones(260,1);

for i = 1:260 
    if strcmp(ii(i), 'i')
        labels(i) = 1;
    elseif strcmp(ii(i), 'ii')
        labels(i) = 2;
    elseif strcmp(ii(i), 'iii')
        labels(i) = 3;
    else
        labels(i) = 4;
    end        
end

[rank_fpqm_scr, rank_fpqm_cho, rank_fpqm_scr_supl] = rank_scr_cho_supl(fpqm, labels);
[rank_fpqm_log_scr, rank_fpqm_log_cho, rank_fpqm_log_scr_supl] = rank_scr_cho_supl(fpqm_log, labels);
[rank_nt_scr, rank_nt_cho, rank_nt_scr_supl] = rank_scr_cho_supl(nt, labels);
[rank_vs_scr, rank_vs_cho, rank_vs_scr_supl] = rank_scr_cho_supl(vs, labels);

csvwrite('rank_fpqm_2.csv', rank_fpqm_scr);
csvwrite('rank_fpqm_cho.csv', rank_fpqm_cho);
csvwrite('rank_fpqm_1.csv', rank_fpqm_scr_supl);

csvwrite('rank_fpqm_log_2.csv', rank_fpqm_log_scr);
csvwrite('rank_fpqm_log_cho.csv', rank_fpqm_log_cho);
csvwrite('rank_fpqm_log_1.csv', rank_fpqm_log_scr_supl);

csvwrite('rank_vs_2.csv', rank_vs_scr);
csvwrite('rank_vs_cho.csv', rank_vs_cho);
csvwrite('rank_vs_1.csv', rank_vs_scr_supl);

csvwrite('rank_nt_2.csv', rank_nt_scr);
csvwrite('rank_nt_cho.csv', rank_nt_cho);
csvwrite('rank_nt_1.csv', rank_nt_scr_supl);
