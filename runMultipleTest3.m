function runMultipleTest3(result_file)
    %pos_files={'EWSFLI_ews502_hg18_10000','EWSFLI_Huvec_hg18_10000','GR_3134_mm8_10000','GR_att20_mm8_10000','ESRRB_mm8_10000'};
    %pos_files={'ESRRB_mm8_10000','GR_3134_mm8_10000','GR_att20_mm8_10000'};
    %pos_files={'2_HNF4A_prolCell_pos','3_GATA6_prolCell_pos','4_CDX2_prolCell_pos'};
    %pos_files={'1_TAL1_erythroid_pos','10_H3K4me2_prolCell_pos'};
    pos_files={'2_HNF4A_prolCell_pos','3_GATA6_prolCell_pos','4_CDX2_prolCell_pos','1_TAL1_erythroid_pos','10_H3K4me2_prolCell_pos'};
    %pos_files={'2_HNF4A_diffCell_pos','3_GATA6_diffCell_pos','4_CDX2_diffCell_pos','1_TAL1_jurkat_pos','10_H3K4me2_diffCell_pos'};
    %pos_files={'2_HNF4A_prolCell_pos'};
    %pos_files={'4_CDX2_prolCell_pos','4_CDX2_diffCell_pos','1_TAL1_erythroid_pos','1_TAL1_jurkat_pos'};
    %pos_files={'ESRRB_mm8_10000','GR_3134_mm8_10000','GR_att20_mm8_10000'};
    %pos_files={'EWSFLI_Huvec_hg18_10000','EWSFLI_ews502_hg18_10000'};
    %pos_files={'EWSFLI_Huvec_hg18_10000'};
    %pos_files={'test1'};
    %neg_files={'EWSFLI_ews502_hg18_neg10x_10000','EWSFLI_Huvec_hg18_neg10x_10000','GR_3134_mm8_neg10x_10000','GR_att20_neg10x_10000','ESRRB_mm8_neg10x_10000'};
    %neg_files={'ESRRB_mm8_neg10x_10000','GR_3134_mm8_neg10x_10000','GR_att20_neg10x_10000'};
    %neg_files={'ESRRB_mm8_neg10x_10000','GR_3134_mm8_neg10x_10000','GR_att20_neg10x_10000'};
    %neg_files={'2_HNF4A_prolCell_nullx10','3_GATA6_prolCell_nullx10','4_CDX2_prolCell_nullx10'};
    %neg_files={'1_TAL1_erythroid_nullx10','10_H3K4me2_prolCell_nullx10'};
    neg_files={'2_HNF4A_prolCell_nullx10','3_GATA6_prolCell_nullx10','4_CDX2_prolCell_nullx10','1_TAL1_erythroid_nullx10','10_H3K4me2_prolCell_nullx10'};
    %neg_files={'2_HNF4A_diffCell_nullx10','3_GATA6_diffCell_nullx10','4_CDX2_diffCell_nullx10','1_TAL1_jurkat_nullx10','10_H3K4me2_diffCell_nullx10'};
    %neg_files={'2_HNF4A_prolCell_nullx10'};
    %neg_files={'4_CDX2_diffCell_pos','4_CDX2_prolCell_pos','1_TAL1_jurkat_pos','1_TAL1_erythroid_pos'};
    %neg_files={'ESRRB_mm8_neg10x_10000','GR_3134_mm8_neg10x_10000','GR_att20_neg10x_10000'};
    %neg_files={'EWSFLI_Huvec_hg18_neg10x_10000','EWSFLI_ews502_hg18_neg10x_10000'};
    %neg_files={'EWSFLI_Huvec_hg18_neg10x_10000'};
    %neg_files={'test2'};
    
    %nos_selected_kmers = [10,50];
    %nos_selected_kmers = [100,50];
    %nos_selected_kmers = [100,50,30,10];
    nos_selected_kmers = [100];
    %nos_selected_kmers = [200];
    %nos_selected_kmers = [100];
    %nos_selected_kmers = [2];
    kmers_length = 6;
    %kmers_length = 4;
    %kmers_length = 2;
    %combination_distances = [10,30,50,100,200];
    combination_distances = [100];
    %combination_distances = [2];
    %svm_kmer_select_methods = {'+','+-','-'};
    %svm_kmer_select_methods = {'3+','3+-','3-'};
    %svm_kmer_select_methods = {'+'};
    %svm_kmer_select_methods = {'2+d'};
    %svm_kmer_select_methods = {'2+o'};
    %svm_kmer_select_methods = {'2+n'};
    svm_kmer_select_methods = {'2+n','2+-n','2-n'};
    %svm_kmer_select_methods = {'2+o','2+-o','2-o'};
    %svm_kmer_select_methods = {'2+d'};
    n_folds = 5;
    
    preselects = {'mRMR','tTest','infoGain','Fisher','adaboost'};
    
    %load done test
    done_tests = load_done_test(result_file);
    
    for fi=1:length(pos_files)
        %Get the data file
        pos_file = pos_files{fi};
        neg_file = neg_files{fi};
        
        %Run test
        %Run benchmark (Lee method)
        %{
        if ~check_if_done(done_tests,'enhancer_svm', pos_file, neg_file,0,kmers_length,0, 'none','none',n_folds)
            disp(['Running: ','enhancer_svm(',pos_file,',',neg_file,',',num2str(kmers_length),',',num2str(n_folds),')']);
            a = enhancer_svm(pos_file,neg_file,kmers_length,n_folds);

            [ROC_AUC,PR_AUC,weight] = a.plot_roc(true);
            disp(['ROC AUC: ',num2str(ROC_AUC)]);
            disp(['PR AUC: ',num2str(PR_AUC)]);

            clear a;

            store_data(result_file,'enhancer_svm', pos_file, neg_file,0,kmers_length,0, 'none','none',n_folds, ROC_AUC,PR_AUC,weight);
        end
        %}
        
        
        %Run new method with different parameter
        for no=1:length(nos_selected_kmers)
            %no slected kmers
            no_selected_kmers = nos_selected_kmers(no);
            for cd=1:length(combination_distances)
                %Combination distance
                combination_distance = combination_distances(cd);
                
                %{
                for prs=1:length(preselects)
                    %preselect other
                    preselect = preselects{prs};

                    if ~check_if_done(done_tests,'enhancer_combined_svm', pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,'none',n_folds)
                        disp(['Running: ','enhancer_combined_svm(',pos_file,',',neg_file,',',num2str(no_selected_kmers),',',num2str(kmers_length),',',num2str(combination_distance),',',preselect,',n,',num2str(n_folds),')']);
                        a = enhancer_combined_svm(pos_file,neg_file,no_selected_kmers,kmers_length,combination_distance,preselect,'n',n_folds);

                        [ROC_AUC,PR_AUC,weight] = a.plot_roc(true);
                        disp(['ROC AUC: ',num2str(ROC_AUC)]);
                        disp(['PR AUC: ',num2str(PR_AUC)]);

                        clear a;

                        store_data(result_file,'enhancer_combined_svm', pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,'n',n_folds, ROC_AUC,PR_AUC,weight);
                    end
                end
                %}
                

                %preselect svm
                preselect='svm';
                
                %{}
                for sm=1:length(svm_kmer_select_methods)
                    %svm preselect method
                    svm_kmer_select_method = svm_kmer_select_methods{sm};
                    
                    if ~check_if_done(done_tests,'enhancer_combined_svm', pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,svm_kmer_select_method,n_folds)
                        disp(['Running: ','enhancer_combined_svm(',pos_file,',',neg_file,',',num2str(no_selected_kmers),',',num2str(kmers_length),',',num2str(combination_distance),',',preselect,',',svm_kmer_select_method,',',num2str(n_folds),')']);
                        a = enhancer_combined_svm(pos_file,neg_file,no_selected_kmers,kmers_length,combination_distance,preselect,svm_kmer_select_method,n_folds);

                       [ROC_AUC,PR_AUC,weight] = a.plot_roc(true);
                        disp(['ROC AUC: ',num2str(ROC_AUC)]);
                        disp(['PR AUC: ',num2str(PR_AUC)]);

                        clear a;

                        store_data(result_file,'enhancer_combined_svm', pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,svm_kmer_select_method,n_folds, ROC_AUC,PR_AUC,weight);
                    end
                    
                    %{
                    if ~check_if_done(done_tests,'enhancer_combined_only_svm', pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,svm_kmer_select_method,n_folds)
                        disp(['Running: ','enhancer_combined_only_svm(',pos_file,',',neg_file,',',num2str(no_selected_kmers),',',num2str(kmers_length),',',num2str(combination_distance),',',preselect,',',svm_kmer_select_method,',',num2str(n_folds),')']);
                        a = enhancer_combined_only_svm(pos_file,neg_file,no_selected_kmers,kmers_length,combination_distance,preselect,svm_kmer_select_method,n_folds);

                       [ROC_AUC,PR_AUC,weight] = a.plot_roc(true);
                        disp(['ROC AUC: ',num2str(ROC_AUC)]);
                        disp(['PR AUC: ',num2str(PR_AUC)]);

                        clear a;

                        store_data(result_file,'enhancer_combined_only_svm', pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,svm_kmer_select_method,n_folds, ROC_AUC,PR_AUC,weight);
                    end
                    %}
                    
                    %{
                    if ~check_if_done(done_tests,'enhancer_combined_without_top_svm', pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,svm_kmer_select_method,n_folds)
                        disp(['Running: ','enhancer_combined_without_top_svm(',pos_file,',',neg_file,',',num2str(no_selected_kmers),',',num2str(kmers_length),',',num2str(combination_distance),',',preselect,',',svm_kmer_select_method,',',num2str(n_folds),')']);
                        a = enhancer_combined_only_svm(pos_file,neg_file,no_selected_kmers,kmers_length,combination_distance,preselect,svm_kmer_select_method,n_folds);

                       [ROC_AUC,PR_AUC,weight] = a.plot_roc(true);
                        disp(['ROC AUC: ',num2str(ROC_AUC)]);
                        disp(['PR AUC: ',num2str(PR_AUC)]);

                        clear a;

                        store_data(result_file,'enhancer_combined_without_top_svm', pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect,svm_kmer_select_method,n_folds, ROC_AUC,PR_AUC,weight);
                    end
                    %}
                end
                %}
                
            end
        end
    end
end

function done_tests = load_done_test(file)
    %file = 'Test_Result.csv';
    done_tests = {};
    if exist(file,'file') == 2
        fid = fopen(file);
        tline = fgetl(fid);
        tline = fgetl(fid);
        i = 1;
        while ischar(tline)
            done_tests{i} = tline;
            
            i = i + 1;
            tline = fgetl(fid);
        end
        fclose(fid);
    end
end

function done = check_if_done(done_tests,method, pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold)
    for i=1:length(done_tests)
        done_test = done_tests{i};
        test = [method,',',pos_file,',',neg_file,',',num2str(no_selected_kmers),',',num2str(kmers_length),',',num2str(combination_distance),',',preselect_method,',',svm_select_method,',',num2str(no_fold)];
        
        if(~isempty(strfind(done_test,test)))
            done = true;
            return;
        end
    end
    done = false;
end

function store_data(file,method, pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold, ROC_AUC, PR_AUC,weight)
    %file = 'Test_Result.csv';
    %Write result
    if exist(file,'file') == 2
        fid = fopen(file,'a');
        fprintf(fid,'%s,%s,%s,%d,%d,%d,%s,%s,%d,%f,%f\n',method, pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold, ROC_AUC,PR_AUC);
        fclose(fid);
    else
        fid = fopen(file,'w');
        fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n','method', 'pos_file', 'neg_file','no_selected_kmers','kmers_length','combination_distance',' preselect_method','svm_select_method','no_fold','ROC AUC','PR AUC');
        fprintf(fid,'%s,%s,%s,%d,%d,%d,%s,%s,%d,%f,%f\n',method, pos_file, neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold, ROC_AUC,PR_AUC);
        fclose(fid);
    end
    
    %Store weight
    weight_file = ['weights/',method,'_',pos_file,'_',neg_file,'_',num2str(kmers_length),'_combined_',num2str(combination_distance),'_',preselect_method,svm_select_method,'_selected_',num2str(no_selected_kmers),'.tsv'];
    store_weight(weight_file,weight)
end

function store_weight(file,weight)
    fid = fopen(file,'w');
    for i=1:size(weight,1);
        fprintf(fid,'%s\t%f\n',weight{i,1},weight{i,2});
    end
    fclose(fid);
end
