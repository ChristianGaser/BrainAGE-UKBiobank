function [num_array, age, male, BA, BA_weighted] = cg_analyze_ukb_long(codes, analysis, n_covariates, add_nuisance, order_poly, names, BA_weighting, force_regres, data_dir)
% Prepare analysis of longitudinal UKB data and call VBM or BA analysis
% Always call that function from the analysis folder and organize your data and folders in the following way:
%  your_data_folder/analysis
%  your_data_folder/tables
%  your_data_folder/diff_s6mwmwp1r
%  your_data_folder/s6mwmwp1r
%  your_data_folder/s6mwmwp2r
%
% Format [num_array, age, male, BA, BA_weighted] = cg_analyze_ukb_long(codes, analysis, n_covariates, add_nuisance, order_poly, names, BA_weighting, force_regres, data_dir)
%
% Defined codes are searched in the xls-file which is located in the 'tables' folder in data_dir.
%
% Input:
% codes        - Codes separated by spaces (start with Covariates of interest)
% analysis     - type of analysis: 0 - Check UKB and return variables
%                                  1 - VBM GM Baseline
%                                  2 - VBM WM Baseline
%                                  3 - VBM GM difference Follow up - Baseline
%                                  4 - BrainAGE Baseline
%                                  5 - BrainAGE diff Follow up - Baseline
%                                  6 - BrainAGE both
% n_covariates - number of covariates of interest to build contrasts for VBM
% add_nuisance - add nuisance:     0 - None
%                                  1 - Add age
%                                  2 - Add sex
%                                  3 - Add age+sex
% order_poly   - order of polynomial fit for BA analysis
% names        - char array of variable names
% BA_weighting - use median BA or other weightings instead of GLM-weighting for all 8 models
%                  0 - use all 8 models and estimate weighting by maximizing variance to defined parameter
%                  1 - use median of 8 models
%                  2 - use GLM estimation to estimate model weights to minimize MAE
%                  3 - use RVR to combine models 
% force_regres - force regression (helpful if categorical variables are not fully recognized)
% data_dir     - data directory with data subfolders and analysis subfolder
%
% Output:
% The returned parameters are limited to those subjects where all codes were correctly
% defined (i.e. no NANs, no negative values for categorical parameters)
%
% num_array   - returned number array for selected codes and subjects
% age         - age of selected subjects
% sex         - sex of selected subjects (0 - female, 1 - male)
% BA          - BA for all 8 models of selected subjects (or 4 weighted models for BA_weighted > 0)
% BA_weighted - GLM-weighted BA of selected subjects (variance maximized for regression)

% data subfolders
diff_mwmwp1r_dir = 'diff_mwmwp1r_r1900';
mwmwp1r_dir      = 'mwmwp1r_r1900_2';
mwmwp2r_dir      = 'mwmwp2r_r1900_2';
tables_dir       = 'tables';

% filenames for data
diff_mwmwp1r = 'diff_s6mwmwp1r';
mwmwp1r      = 's6mwmwp1r';
mwmwp2r      = 's6mwmwp2r';
age2_file    = 'age-2.txt';
age3_file    = 'age-3.txt';
male_file    = 'male.txt';
TIV2_file    = 'TIV-2.txt';
GM_mask      = 'GM_015.nii';
WM_mask      = 'WM_015.nii';

use_GLM = 1; % alternative is optimization, with forcing pos. weights

BA2_file     = 'BA-2_8models_GPR.txt';
BA3_file     = 'BA-3_8models_GPR.txt';
BAW2_file    = 'BA-2_5weightings.txt';
BAW3_file    = 'BA-3_5weightings.txt';

% outputs
BA           = [];
BA_weighted  = [];

if nargin < 1
  codes = spm_input('Codes separated by spaces (start with Covs):',1,'s');
end

if nargin < 2
  analysis = spm_input('Which analysis?','+1','m','Check UKB and return variables|GM Baseline|WM Baseline|GMdiff Follow up - Baseline|BrainAGE Baseline|BrainAGE diff Follow up - Baseline|BA both',0:6,2);
end

% codes can be also a matrix and is then straightforward used as number array
% without loading the excel-sheet
if isnumeric(codes) || islogical(codes)
  num_array  = double(codes);
  num_array1 = double(codes(:,1));
  [x,y] = size(num_array);
  [n_parameters, ind] = min([x y]);
  
  % check whether array is continuous
  n_categories = numel(unique(num_array1(isfinite(num_array1) & num_array1 >= 0)));
  for i=1:size(num_array,2)
    num_arrayi = num_array(:,i);
    if sum(num_arrayi(isfinite(num_arrayi)) - round(num_arrayi(isfinite(num_arrayi)))) == 0 && n_categories < 10
      is_continuous{i} = 0;
    else
      is_continuous{i} = 1;
    end  
  end
  
  % check whether there are negative categorical numbers that should be excluded
  if min(num_array(isfinite(num_array))) < 0 && numel(unique(num_array1(num_array1 < 0))) < 5
    % remove negative values for categorical data and keep finite data
    if n_parameters == 1
      finite_rows = isfinite(num_array) & num_array >= 0;
    else
      finite_rows = all(isfinite(num_array')) & all(num_array' >= 0);
    end
  else
    if n_parameters == 1
      finite_rows = isfinite(num_array);
    else
      finite_rows = all(isfinite(num_array'));
    end
  end  

  if ind == 1
    num_array = num_array';
  end

  ismatrix_codes = 1;
  if nargin < 6
    names = strsplit(num2str(1:n_parameters),' ');
  end
else
  ismatrix_codes = 0;
end

% convert codes if defined as characters with spaces as delimiter
if ~ismatrix_codes && ischar(codes)
  codes = strsplit(codes,' ');
end

if ~ismatrix_codes, n_parameters = numel(codes); end

% define number of covariates for analyis
if ~ismatrix_codes && nargin < 3
  if n_parameters > 1 && analysis
    n_covariates = spm_input('How many covariates of interest:','+1','n',1);
  else
    n_covariates = 1;
  end
end

if nargin < 4
  add_nuisance = spm_input('Add nuisance?',1,'m','None|Add age|Add sex|Add age+sex',0:3,0);
end

if nargin < 5
  if analysis > 3
    order_poly = spm_input('Polynomial order?','+1','n',1);
  else
    order_poly = 1;
  end
end

% use median BA or all 8 models
if nargin < 7
  BA_weighting = 0;
end

if nargin < 8
  force_regres = 0;
end

% either define here the full path where your data are located or use relative path
if nargin < 9
  data_dir         = '../';
end

% check whether all folders exist
if analysis && analysis < 4
  if ~exist(fullfile(data_dir,tables_dir),'dir')
    if exist(fullfile('.',tables_dir),'dir')
      data_dir = '.';
    else
      error('Table folder %s does not exist. Folder structure should be:\nyour_data_folder/analysis\nyour_data_folder/tables\nyour_data_folder/diff_s6mwmwp1r\nyour_data_folder/s6mwmwp1r\nyour_data_folder/s6mwmwp2r\n',fullfile(data_dir,tables_dir));
    end
  end
  
  if ~exist(fullfile(data_dir,tables_dir,'ukb49261_long3046.xlsx'),'file')
    error('Excel file %s does not exist. Folder structure should be:\nyour_data_folder/analysis\nyour_data_folder/tables\nyour_data_folder/diff_s6mwmwp1r\nyour_data_folder/s6mwmwp1r\nyour_data_folder/s6mwmwp2r\n',fullfile(data_dir,tables_dir,'ukb49261_long3046.xlsx'));
  end
  
  if ~exist(fullfile(data_dir,diff_mwmwp1r_dir),'dir')
    error('Data folder %s does not exist. Folder structure should be:\nyour_data_folder/analysis\nyour_data_folder/tables\nyour_data_folder/diff_s6mwmwp1r\nyour_data_folder/s6mwmwp1r\nyour_data_folder/s6mwmwp2r\n',fullfile(data_dir,diff_mwmwp1r));
  end
end

if ~ismatrix_codes
  [num_array, names, finite_rows] = cg_get_ukb_data(fullfile(data_dir,tables_dir,'ukb49261_long3046.xlsx'), codes);

  % check whether array is continuous
  n_categories = numel(unique(num_array(isfinite(num_array) & num_array >= 0)));
  for i=1:size(num_array,2)
    num_arrayi = num_array(:,i);
    if sum(num_arrayi(isfinite(num_arrayi)) - round(num_arrayi(isfinite(num_arrayi)))) == 0 && n_categories < 10
      is_continuous{i} = 0;
    else
      is_continuous{i} = 1;
    end  
  end 
end

% get age at BL and FU
TIV2 = spm_load(fullfile(data_dir,tables_dir,TIV2_file));
age2 = spm_load(fullfile(data_dir,tables_dir,age2_file));
age3 = spm_load(fullfile(data_dir,tables_dir,age3_file));
male = spm_load(fullfile(data_dir,tables_dir,male_file));

if BA_weighting
  BA2_name = BAW2_file;
  BA3_name = BAW3_file;
else
  BA2_name = BA2_file;
  BA3_name = BA3_file;
end

% build parameters w.r.t. selected analyis
switch analysis
  case 1, subdir = mwmwp1r_dir;      pattern = mwmwp1r;      mask = GM_mask; gsc = TIV2; gsc_str = 'TIV'; % use GM mask and TIV for global scaling
  case 2, subdir = mwmwp2r_dir;      pattern = mwmwp2r;      mask = WM_mask; gsc = TIV2; gsc_str = 'TIV';  % use WM mask and TIV for global scaling
  case 3, subdir = diff_mwmwp1r_dir; pattern = diff_mwmwp1r; mask = GM_mask; gsc = age3-age2; gsc_str = 'age';  % use GM mask and age-diff for global scaling
  case 4, P = spm_load(fullfile(data_dir,tables_dir,BA2_name));
  case 5, P2 = spm_load(fullfile(data_dir,tables_dir,BA2_name)); P3 = spm_load(fullfile(data_dir,tables_dir,BA3_name)); P = P3 - P2;
  case 6, P2 = spm_load(fullfile(data_dir,tables_dir,BA2_name)); P3 = spm_load(fullfile(data_dir,tables_dir,BA3_name));
end

% find files
if analysis < 4 && analysis
  if ismatrix_codes
    str_name = names;
  else
    str_name = codes;
  end
  analysis_dir = [pattern '_' str_name{1}];
  for i=2:numel(str_name)
    analysis_dir = [analysis_dir '_' str_name{i}];
  end
  P = spm_select('FPList',fullfile(data_dir,subdir),pattern);
end

% only use rows where all numbers are finite
fprintf('Use %d with finite and positive numbers\n',sum(finite_rows));
num_array = num_array(finite_rows,:);
age2 = age2(finite_rows);
age3 = age3(finite_rows);
male = male(finite_rows);

if analysis
  if analysis < 6
    P = P(finite_rows,:);
    if analysis == 5, P2 = P2(finite_rows,:); end
  else
    P = [P2(finite_rows,:); P3(finite_rows,:)];
    num_array = [num_array; num_array];
  end
  if analysis < 4
    gsc = gsc(finite_rows); 
  else
    % rescue BrainAge
    if analysis > 3
      BA = P;
      if analysis == 5, BA2 = P2; end
    end
  end
end

% add predefined nuisance parameters
switch add_nuisance
  case 1, num_array = [num_array age2];      n_parameters = n_parameters + 1; names{end+1} = 'age'; is_continuous{end+1} = 1;
  case 2, num_array = [num_array male];      n_parameters = n_parameters + 1; names{end+1} = 'male'; is_continuous{end+1} = 0;
  case 3, num_array = [num_array age2 male]; n_parameters = n_parameters + 2; names{end+1} = 'age'; names{end+1} = 'male'; is_continuous{end+1} = 1; is_continuous{end+1} = 0;
end

% if number of categories is between 4..12 we have to ask how to deal with that
if ~is_continuous{1}
  n_categories = numel(unique(num_array(:,1)));
  if n_categories > 3 && ~force_regres
    figure(11)
    cat_plot_boxplot(num_array(:,1),struct('style',0,'showdata',2));
%    title(names{1})
    is_continuous{1} = spm_input('Should we apply regression or Anova for group differences?',1,'m','Anova|Regression',0:1,0);
  elseif ~force_regres 
    is_continuous{1} = 0;
  else
    is_continuous{1} = 1; % only force regression for 1st parameter
  end
end

if n_parameters > 1
  cc = corrcoef(num_array);
  fprintf('Correlations between your covariates of interest and nuisance parameters:\n');
  disp(triu(cc))
end

% use multiple regression with global scaling (TIV or age difference) and estimate model
if analysis < 4 && analysis % VBM analysis

  fprintf('Run analysis of %s with %s as global scaling and %d covariates\n',pattern,gsc_str,n_parameters)
  matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(fullfile(data_dir,'analysis',analysis_dir));
  matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = cellstr(P);
  matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});
  matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
  
  % define covariates
  for i=1:n_parameters
    matlabbatch{1}.spm.stats.factorial_design.cov(i).c = num_array(:,i);
    matlabbatch{1}.spm.stats.factorial_design.cov(i).cname = names{i};
    matlabbatch{1}.spm.stats.factorial_design.cov(i).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(i).iCC = 5;
  end
  
  matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
  matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
  matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
  matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(mask);
  matlabbatch{1}.spm.stats.factorial_design.globalc.g_user.global_uval = gsc;
  matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_yes.gmscv = mean(gsc);
  matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 2;
  
  % estimate model
  matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
  matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
  matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
  
  % define contrasts
  matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
  matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'effects of interest';
  matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [zeros(n_covariates) eye(n_covariates)];
  matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
  matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'positive regression';
  matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1];
  matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
  matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'negative regression';
  matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 -1];
  matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
  matlabbatch{3}.spm.stats.con.delete = 1;  
  
  % save batch for later edits
  save analysis_dir matlabbatch
  
  % run job
  spm_jobman('run',matlabbatch);
end

% show boxplots
if ~analysis
  for i = n_parameters:-1:1
    figure(10+i)
    cat_plot_boxplot(num_array(:,i),struct('style',0,'showdata',2))
    if i > 1, title(names{i}); end
  end
elseif analysis > 3 % BA-analysis
  for i = n_parameters:-1:1
    % make multiple regression
    if is_continuous{i}
      X = cat_stat_polynomial(num_array(:,i),order_poly);
      switch BA_weighting
        case 0
          Y = BA;
          if use_GLM % use GLM
            Beta = pinv(Y)*X;
            Beta = Beta./sum(Beta);
            BA_weighted = sum(BA*Beta,2);
    
            % scale estimated BA values by ratio between SD of original and estimated BA values to get the same range
            BA_weighted = BA_weighted*mean(std(BA))/mean(std(BA_weighted));
            
            if analysis == 5
              BA_weighted2 = sum(BA2*Beta,2);
              BA_weighted2 = BA_weighted2*mean(std(BA2))/mean(std(BA_weighted2));
            end
          else
            num_ensembles = size(Y,2);
        
            % Objective function: Mean Squared Error of weighted predictions
            objective = @(weights) mean((X - Y * weights).^2);
            
            % Initial guess
            initial_weights = ones(num_ensembles, 1) / num_ensembles;
        
            % Linear equality constraint to make weights sum to 1
            Aeq = ones(1, num_ensembles);
            beq = 1;
        
            % Bounds to ensure weights are non-negative
            lb = zeros(num_ensembles, 1) + 0.001;
            ub = ones(num_ensembles, 1) - 0.001;
        
            % Options: Display iterations
            options = optimoptions('fmincon', 'Display', 'iter');
        
            % Optimization
            [optimal_weights, fval] = fmincon(objective, initial_weights, [], [], Aeq, beq, lb, ub, [], options)
            BA_weighted = sum(BA*optimal_weights,2);

          end
          
          Y = BA_weighted;
        case 1
          % use median of all BA models
          Y = BA(:,3);
        case 2
          % use GLM to minimize MAE
          Y = BA(:,2);
        case 3
          % use RVR to minimize RVR
          Y = BA(:,4);
        case 4
          % use RVR to minimize RVR
          Y = BA(:,5);
      end
      
      % remove effects due to additional nuisance paramaters
      if size(num_array,2) - n_covariates
        % also add baseline BA to nuisance for analysis 5
        if (add_nuisance == 4) && (analysis == 5)
          G = [ones(size(num_array,1),1) num_array(:,n_covariates+1:end) BA_weighted2];
        else
          G = [ones(size(num_array,1),1) num_array(:,n_covariates+1:end)];
        end
        Beta = pinv(G)*Y;
        
        if exist('fitlm') && i==1
          if (add_nuisance == 4) % remove last entry which is BL-BA
            mdl = fitlm([num_array(:,i) G(:,2:end-1)],Y)
          else
            mdl = fitlm([num_array(:,i) G(:,2:end)],Y)
          end
        end
        
        % and remove effects for all data
        Y = Y - G*Beta;
      end

      % correct BrainAGE difference w.r.t. time gap
      if analysis == 5
        Y = Y./(age3-age2) * mean(age3-age2);
        fprintf('Correct BrainAGE difference w.r.t. time gap of %g years (SD=%g).\n',mean(age3-age2),std(age3-age2));
      end
      
      fprintf('%s\n',names{i});
      cat_plot_scatter(num_array(:,i),Y,'fig',10+i,'fit_poly',order_poly,'jitter',1)
      fprintf('\n')
      fprintf('Mean/Median (Std) for parameter %s (n = %d): %g/%g (%g)\n',names{i},numel(num_array(:,i)),mean(num_array(:,i)),median(num_array(:,i)),std(num_array(:,i)));
      if i > 1, title(names{i}); end
    else
      
      if n_covariates > 1
        error('For Anova with groups only one variable is allowed.');
      end
      categories = unique(num_array(:,i));

      c = eye(numel(categories));
      contrast = [];
      ind = [];
      for j=1:numel(categories)
        ind_groups{j} = find(num_array(:,i) == categories(j));
        contrast = [contrast; repmat(c(j,:),numel(ind_groups{j}),1)];
        ind = [ind; ind_groups{j}];
      end
      % re-order w.r.t. ind_groups
      if max(ind) == numel(ind)
        contrast = contrast(ind,:);
      else
        fprintf('Warning: Order of ind_groups should be not mixed because we cannot identify the correct order!\n');
      end

      X = contrast;
      if analysis == 4
        % use median of all BA models
        Y = median(BA,2);
      else
        Y = BA;
        Beta = pinv(Y)*X;
        Beta = Beta./sum(Beta);
        BA_weighted = sum(BA*Beta,2);

        % scale estimated BA values by ratio between SD of original and estimated BA values to get the same range
        BA_weighted = BA_weighted*mean(std(BA))/mean(std(BA_weighted));
        Y = BA_weighted;
      end

      % remove effects due to additional nuisance paramaters
      if size(num_array,2) - n_covariates
        G = [ones(size(num_array,1),1) num_array(:,n_covariates+1:end)];

        Beta = pinv(G)*Y;

        % and remove effects for all data
        Y = Y - G*Beta;
      end

      % correct BrainAGE difference w.r.t. time gap
      if analysis == 5
        Y = Y./(age3-age2) * mean(age3-age2); 
        fprintf('Correct BrainAGE difference w.r.t. time gap of %g years (SD=%g).\n',mean(age3-age2),std(age3-age2));
      end

      Ytmp = cell(numel(categories),1);
      for j = 1:numel(categories)
        num_arrayi = num_array(:,i);
        num_arrayi = num_arrayi(isfinite(num_arrayi));
        Ytmp{j} = Y(num_arrayi == categories(j));
        fprintf('Mean/Median (SD) for category %d (n = %d): %g/%g (%g)\n',categories(j),numel(Ytmp{j}),mean(Ytmp{j}),median(Ytmp{j}),std(Ytmp{j}));
      end
      
      try
        figure(10+i)
      cat_plot_boxplot(Ytmp,struct('style',0,'showdata',2,'names',num2str(categories)));
      if i > 1, title(names{i}); end
      end
    end
    
  end

  BA_weighted = Y;
end

age = age2;
fprintf('Age BL: %g years (SD=%g).\n',mean(age2),std(age2));
fprintf('Age FU: %g years (SD=%g).\n',mean(age3),std(age3));

if ~nargout
  clear num_array age male BA BA_weighted
end