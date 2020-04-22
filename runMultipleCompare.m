function runMultipleTest()
    %pos_files={'EWSFLI_ews502_hg18_10000','EWSFLI_Huvec_hg18_10000','GR_3134_mm8_10000','GR_att20_mm8_10000','ESRRB_mm8_10000'};
    pos_files={'EWSFLI_Huvec_hg18_10000','EWSFLI_ews502_hg18_10000'};
    %pos_files={'EWSFLI_Huvec_hg18_10000'};
    %pos_files={'test1'};
    %pos_files={'test3'};
    %neg_files={'EWSFLI_ews502_hg18_neg10x_10000','EWSFLI_Huvec_hg18_neg10x_10000','GR_3134_mm8_neg10x_10000','GR_att20_neg10x_10000','ESRRB_mm8_neg10x_10000'};
    neg_files={'EWSFLI_Huvec_hg18_neg10x_10000','EWSFLI_ews502_hg18_neg10x_10000'};
    %neg_files={'EWSFLI_Huvec_hg18_neg10x_10000'};
    %neg_files={'test2'};
    %neg_files={'test4'};
    
    %nos_selected_kmers = [10,30,50,100];
    %nos_selected_kmers = [200];
    nos_selected_kmers = [100];
    %nos_selected_kmers = [2];
    kmers_length = 6;
    %kmers_length = 4;
    %kmers_length = 2;
    %combination_distances = [10,30,50,100];
    combination_distances = [100];
    %combination_distances = [10,50,100];
    %combination_distances = [2];
    %svm_kmer_select_methods = {'+','+-','-'};
    %svm_kmer_select_methods = {'3+','3+-','3-'};
    %svm_kmer_select_methods = {'+'};
    %svm_kmer_select_methods = {'2+d'};
    svm_kmer_select_methods = {'2+o'};
    %svm_kmer_select_methods = {'2+d'};
    n_folds = 5;
    n_tries = 5;
    
    for fi=1:length(pos_files)
        %Get the data file
        pos_file = pos_files{fi};
        neg_file = neg_files{fi};
        
        %Run new method with different parameter
        for no=1:length(nos_selected_kmers)
            %no slected kmers
            no_selected_kmers = nos_selected_kmers(no);
            for cd=1:length(combination_distances)
                %Combination distance
                combination_distance = combination_distances(cd);
                 
                %preselect svm
                preselect='svm';
                
                for sm=1:length(svm_kmer_select_methods)
                    %svm preselect method
                    svm_kmer_select_method = svm_kmer_select_methods{sm};
                    
                    if ~check_if_done(pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,svm_kmer_select_method,n_folds,n_tries)
                        disp(['Running: ','combined_vs_lee_svm(',pos_file,',',neg_file,',',num2str(no_selected_kmers),',',num2str(kmers_length),',',num2str(combination_distance),',',preselect,',',svm_kmer_select_method,',',num2str(n_folds),',',num2str(n_tries),')']);
                        a = combined_vs_lee_svm(pos_file,neg_file,no_selected_kmers,kmers_length,combination_distance,preselect,svm_kmer_select_method,n_folds,n_tries);

                        disp('ROC AUC:')
                        ROC_AUC = a.ROC_AUC
                        disp('PR_AUC:')
                        PR_AUC = a.PR_AUC

                        clear a;

                        store_data(pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,svm_kmer_select_method,n_folds,n_tries,ROC_AUC,PR_AUC);
                    end
                end
            end
        end
    end
end

function done = check_if_done(pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold,no_tries)
    file = ['compare/',pos_file,'_',neg_file,'_',num2str(no_selected_kmers),'_',num2str(kmers_length),'_',num2str(combination_distance),'_',preselect_method,'_',svm_select_method,'_',num2str(no_fold),'_',num2str(no_tries)];
    done = exist(file,'file') == 2;
end

function store_data(pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold,no_tries, ROC_AUC, PR_AUC)
    file = ['compare/',pos_file,'_',neg_file,'_',num2str(no_selected_kmers),'_',num2str(kmers_length),'_',num2str(combination_distance),'_',preselect_method,'_',svm_select_method,'_',num2str(no_fold),'_',num2str(no_tries)];
    
    %Write result
    fid = fopen(file,'w');
    fprintf(fid,['combined_vs_lee_svm(',pos_file,',',neg_file,',',num2str(no_selected_kmers),',',num2str(kmers_length),',',num2str(combination_distance),',',preselect_method,',',svm_select_method,',',num2str(no_fold),',',num2str(no_tries),')']);
    
    fprintf(fid,'\n\nROC_AUC:\n');
    for i = 1:size(ROC_AUC,1)
        fprintf(fid,'%f\t%f\n',ROC_AUC(i,1),ROC_AUC(i,2));
    end

    fprintf(fid,'\n\nPR_AUC:\n');
    for i = 1:size(PR_AUC,1)
        fprintf(fid,'%f\t%f\n',PR_AUC(i,1),PR_AUC(i,2));
    end
end
