function summary(file,ofile)
    data = readData(file);
    
    table = preselectTable(data,100,100);
    exportData(ofile,'Preselection N-100 D-100',table);
    
    table = NDTable(data,'svm','2+n');
    exportData(ofile,'ND svm 2+n',table);
    
    table = NDTable(data,'tTest','2+n');
    exportData(ofile,'ND tTest 2+n',table);
    
    table = BestNDTable(data,'svm','2+n');
    exportData(ofile,'Best ND svm 2+n',table);
    
    table = compareTable(data,100,100,'svm','2+n');
    exportData(ofile,'Compare N-100 D-100 svm 2+n',table);
    
    table = compareTable(data,100,100,'tTest','2+n');
    exportData(ofile,'Compare N-100 D-100 tTest 2+n',table);
end

function [preselect_table] = preselectTable(data,n,d)
    %method,pos_file,neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold,ROC AUC,PR AUC
    
    %preselection method
    gc = {'no_selected_kmers',num2str(n),'combination_distance',num2str(d)};
    rc = {%
        {'pos_file','1_TAL1_erythroid_pos'},
        {'pos_file','1_TAL1_jurkat_pos'},
        {'pos_file','2_HNF4A_prolCell_pos'},
        {'pos_file','2_HNF4A_diffCell_pos'},
        {'pos_file','3_GATA6_prolCell_pos'},
        {'pos_file','3_GATA6_diffCell_pos'},
        {'pos_file','4_CDX2_prolCell_pos'},
        {'pos_file','4_CDX2_diffCell_pos'},
        {'pos_file','10_H3K4me2_prolCell_pos'},
        {'pos_file','10_H3K4me2_diffCell_pos'}
        };
    cc = {%
        {'preselect_method','svm','svm_select_method','2+n'},
        {'preselect_method','svm','svm_select_method','2+-n'},
        {'preselect_method','svm','svm_select_method','2-n'},
        {'preselect_method','adaboost','svm_select_method','2n'},
        {'preselect_method','Fisher','svm_select_method','2+n'},
        {'preselect_method','infoGain','svm_select_method','2+n'},
        {'preselect_method','mRMR','svm_select_method','2n'},
        {'preselect_method','tTest','svm_select_method','2+n'}
        };
    
    r_header = {
        '','SVM+','','SVM+-','','SVM-','','AdaBoost','','Fisher','','InforGain','','mRMR','','TTest','',;
        '','AUROC','AUPRC','AUROC','AUPRC','AUROC','AUPRC','AUROC','AUPRC','AUROC','AUPRC','AUROC','AUPRC','AUROC','AUPRC','AUROC','AUPRC'
        };
    c_header = {
        'TAL1(e)';
        'TAL1(j)';
        'HNF4A(p)';
        'HNF4A(d)';
        'GATA6(p)';
        'GATA6(d)';
        'CDX2(p)';
        'CDX2(d)';
        'H3K4me2(p)';
        'H3K4me2(d)'
        };
    
    table = genTable(data,gc,rc,cc,{'ROC AUC','PR AUC'},0);
    
    preselect_table = [r_header;[c_header,num2cell(table)]];
end

function [nd_table] = NDTable(data,method,selection)
    %method,pos_file,neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold,ROC AUC,PR AUC
    
    %preselection method
    gc = {'preselect_method',method,'svm_select_method',selection};
    rc = {%
        {'pos_file','1_TAL1_erythroid_pos'},
        {'pos_file','1_TAL1_jurkat_pos'},
        {'pos_file','2_HNF4A_prolCell_pos'},
        {'pos_file','2_HNF4A_diffCell_pos'},
        {'pos_file','3_GATA6_prolCell_pos'},
        {'pos_file','3_GATA6_diffCell_pos'},
        {'pos_file','4_CDX2_prolCell_pos'},
        {'pos_file','4_CDX2_diffCell_pos'},
        {'pos_file','10_H3K4me2_prolCell_pos'},
        {'pos_file','10_H3K4me2_diffCell_pos'}
        };
    cc = {%
        {'no_selected_kmers','10','combination_distance','100'},
        {'no_selected_kmers','30','combination_distance','100'},
        {'no_selected_kmers','50','combination_distance','100'},
        {'no_selected_kmers','100','combination_distance','100'},
        {'no_selected_kmers','100','combination_distance','10'},
        {'no_selected_kmers','100','combination_distance','30'},
        {'no_selected_kmers','100','combination_distance','50'},
        {'no_selected_kmers','100','combination_distance','100'},
        {'no_selected_kmers','100','combination_distance','200'}
        };
    
    r_header = {
        '','N','10','30','50','100','100','100','100','100','100';
        '','D','100','100','100','100','10','30','50','100','200'
        };
    c_header = {
        'TAL1(e)','AUROC';'','AUPRC';
        'TAL1(j)','AUROC';'','AUPRC';
        'HNF4A(p)','AUROC';'','AUPRC';
        'HNF4A(d)','AUROC';'','AUPRC';
        'GATA6(p)','AUROC';'','AUPRC';
        'GATA6(d)','AUROC';'','AUPRC';
        'CDX2(p)','AUROC';'','AUPRC';
        'CDX2(d)','AUROC';'','AUPRC';
        'H3K4me2(p)','AUROC';'','AUPRC';
        'H3K4me2(d)','AUROC';'','AUPRC'
        };
    
    table = genTable(data,gc,rc,cc,{'ROC AUC','PR AUC'},1);
    
    nd_table = [r_header;[c_header,num2cell(table)]];
end

function [nd_table] = BestNDTable(data,method,selection)
    %method,pos_file,neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold,ROC AUC,PR AUC
    
    %preselection method
    gc = {'preselect_method',method,'svm_select_method',selection};
    rc = {%
        {'pos_file','1_TAL1_erythroid_pos'},
        {'pos_file','1_TAL1_jurkat_pos'},
        {'pos_file','2_HNF4A_prolCell_pos'},
        {'pos_file','2_HNF4A_diffCell_pos'},
        {'pos_file','3_GATA6_prolCell_pos'},
        {'pos_file','3_GATA6_diffCell_pos'},
        {'pos_file','4_CDX2_prolCell_pos'},
        {'pos_file','4_CDX2_diffCell_pos'},
        {'pos_file','10_H3K4me2_prolCell_pos'},
        {'pos_file','10_H3K4me2_diffCell_pos'}
        };
    cc = {{}};
    
    r_header = {'','ROC','N','D','PR','N','D'};
   c_header = {
        'TAL1(e)';
        'TAL1(j)';
        'HNF4A(p)';
        'HNF4A(d)';
        'GATA6(p)';
        'GATA6(d)';
        'CDX2(p)';
        'CDX2(d)';
        'H3K4me2(p)';
        'H3K4me2(d)'
        };
    
    roctable = genTable(data,gc,rc,cc,{'ROC AUC','no_selected_kmers','combination_distance'},0);
    prtable = genTable(data,gc,rc,cc,{'PR AUC','no_selected_kmers','combination_distance'},0);
    
    nd_table = [r_header;[c_header,num2cell(roctable),num2cell(prtable)]];
end

function [com_table] = compareTable(data,n,d,method,selection)
    %method,pos_file,neg_file,no_selected_kmers,kmers_length,combination_distance, preselect_method,svm_select_method,no_fold,ROC AUC,PR AUC
    
    %preselection method
    gc = {'no_selected_kmers',num2str(n),'combination_distance',num2str(d),'preselect_method',method,'svm_select_method',selection};
    rc = {%
        {'pos_file','1_TAL1_erythroid_pos'},
        {'pos_file','1_TAL1_jurkat_pos'},
        {'pos_file','2_HNF4A_prolCell_pos'},
        {'pos_file','2_HNF4A_diffCell_pos'},
        {'pos_file','3_GATA6_prolCell_pos'},
        {'pos_file','3_GATA6_diffCell_pos'},
        {'pos_file','4_CDX2_prolCell_pos'},
        {'pos_file','4_CDX2_diffCell_pos'},
        {'pos_file','10_H3K4me2_prolCell_pos'},
        {'pos_file','10_H3K4me2_diffCell_pos'}
        };
    cc = {{}};
    
    r_header = {
        '','DSK','';
        '','AUROC','AUPRC'
        };
    c_header = {
        'TAL1(e)';
        'TAL1(j)';
        'HNF4A(p)';
        'HNF4A(d)';
        'GATA6(p)';
        'GATA6(d)';
        'CDX2(p)';
        'CDX2(d)';
        'H3K4me2(p)';
        'H3K4me2(d)'
        };
    
    table = genTable(data,gc,rc,cc,{'ROC AUC','PR AUC'},0);
    
    com_table = [r_header;[c_header,num2cell(table)]];
end

function [table] = genTable(data,gc,rc,cc,sv,v)
    table = [];
    for i=1:length(rc)
        row = [];
        for j=1:length(cc)
            filter = [gc,rc{i},cc{j}];
            r = search(data,filter,sv);
            if(v)
                row = [row,r'];
            else
                row = [row,r];
            end
        end
        table = [table;row];
    end
end

function best = getBest(r,i)
    if(isempty(r))
        best = [0,0];
    else
        r = cellfun(@str2double,r,'un',0);
        nr = cell2mat(r);
        [~,sortedri] =  sort(nr,'descend');
        best = nr(sortedri(1,i),:);
    end
end

function [best] = search(data,condition_paras, lookup_paras)
    data_size = size(data.rows{1});
    
    found_i = 1:data_size(1);
    for i=1:2:length(condition_paras)
        f_i = ismember(data.header,condition_paras{i});
        col = data.rows{f_i};
        found_i = found_i(ismember(col(found_i),condition_paras{i+1}));
    end
    
    results = {};
    for i=1:length(lookup_paras)
        col_i = ismember(data.header,lookup_paras{i});
        col = data.rows{col_i};
        results = [results,col(found_i)];
    end
    
    best = getBest(results,1);
end

function [data] = readData(file)
    fid = fopen(file);
    
    data.header = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s',1, 'Delimiter', ',');
    data.header = [data.header{:}];
    
    data.rows = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s' , 'Delimiter', ',');
    
    fclose(fid);
end

function exportData(file,title,table)
    if exist(file,'file') == 2
        fid = fopen(file,'a');
        fprintf(fid,'\n');
    else
        fid = fopen(file,'w');
    end
    
    fprintf(fid,'%s\n',title);
    
    for i=1:size(table,1)
        for j=1:size(table,2)
            value = table{i,j};
            if(ischar(value))
                fprintf(fid,'%s,',value);
            else
                fprintf(fid,'%f,',value);
            end
        end
        fprintf(fid,'\n');
    end

    fclose(fid);
end
