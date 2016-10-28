function rank_scr_cho_supl

fpqm =  csvread('../data/fpqm .csv', 1,1);
fpqm_log =  csvread('../data/fpqm_log .csv', 1,1);
nt = csvread('../data/nt .csv',1,1);
vs = csvread('../data/vs .csv',1,1);
fpqm_proc = csvread('data/fpqm_proc.csv',1,1);
vs_proc = csvread('data/vs_proc.csv',1,1);

fileid = fopen('stages.txt', 'r');
labels = textscan(fileid, '%s', 'Delimiter', '\n');


labs = ones(260,1);

for i = 1:260 
    if strcmp(labels(i), 'stage i')
        labs(i) = 1;
    elseif strcmp(labels(i), 'stage ii')
        labs(i) = 2;
    elseif strcmp(labels(i), 'stage iii')
        labs(i) = 3;
    else
        labs(i) = 4;
    end        
end

[rank_fpqm_scr, rank_fpqm_cho, rank_fpqm_scr_supl] = rank_scr_cho_supl(fpqm, labels);
[rank_fpqm_log_scr, rank_fpqm_log_cho, rank_fpqm_log_scr_supl] = rank_scr_cho_supl(fpqm_log, labels);
[rank_nt_scr, rank_nt_cho, rank_nt_scr_supl] = rank_scr_cho_supl(nt, labels);
[rank_vs_scr, rank_vs_cho, rank_vs_scr_supl] = rank_scr_cho_supl(vs, labels);

[rank_vs_proc_scr, rank_vs_proc_cho, rank_vs_proc_supl] = rank_scr_cho_supl(vs_proc, labs);
[rank_fpqm_proc_scr, rank_fpqm_proc_cho, rank_fpqm_proc_supl] = rank_scr_cho_supl(fpqm_proc, labs);

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

csvwrite('rank_fpqm_proc_2.csv', rank_fpqm_proc_scr);
csvwrite('rank_fpqm_proc_cho.csv', rank_fpqm_proc_cho);
csvwrite('rank_fpqm_proc_1.csv', rank_fpqm_proc_supl);

csvwrite('rank_vs_proc_2.csv', rank_vs_proc_scr);
csvwrite('rank_vs_proc_cho.csv', rank_vs_proc_cho);
csvwrite('rank_vs_proc_1.csv', rank_vs_proc_supl);