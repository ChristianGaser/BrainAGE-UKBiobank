function [num_array, age, sex, BA, BA_weighted] = cg_analyze_ukb_long_BAonly(sel, num, header, BA_weighting)
% Call cg_analyze_ukb_long and analyze BA with some predefined codes
% Format [num_array, age, sex, BA, BA_weighted] = cg_analyze_ukb_long_BAonly(sel, num, header)
%
% Input:
% sel    - selection of prepared analysis; use negative values for longitudinal BA analysis
% num    - number array of UKB long data without header
% header - header of UKB long data
% BA_weighting - use median BA or other weightings instead of GLM-weighting for all 8 models
%                  0 - use all 8 models and estimate weighting by maximizing variance to defined parameter (default)
%                  1 - use median of 8 models
%                  2 - use GLM estimation to estimate model weights to minimize MAE
%                  3 - use RVR to combine models 
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
%
% Because reading the huge xls-file takes a lot of time you can set the arguments 
% num and headerwith:
%  xls_file = '../tables/ukb49261_long3046.xlsx';
%  [num, header] = cg_get_ukb_data(xls_file);
%
% The returned coefficients are from the fitted polynomial curve.

% some default parameters that are partially overwritten for selected analyses
n_covariates = 1;     % always just one covariate of interest
add_nuisance = 0;     % don't add age or gender as nuisance
order_poly = 1;       % currently only valid for continious BA
force_regression = 1; % force use of regression instead of categorization

% here we can also add BL-BA for the long models as requested by a mean reviewer
%add_nuisance = 4

if nargin < 4
  BA_weighting = 0;        % use full models with GLM maximizing variance to defined parameter
end

if nargin < 3
  xls_file = '../tables/ukb49261_long3046.xlsx';
  [num, header] = cg_get_ukb_data(xls_file);
end

if nargin < 1
  sel = 1; % use negative values for longitudinal BA values
end

% use negative values for longitudinal BA values
if sel > 0
  analysis = 4;
  fprintf('Baseline: ');
else
  analysis = 5;
  fprintf('Longitudinal FU-BL: ');
end

% 3581-2.0 Age at menopause
% 2714-2.0 Age when periods started
% 2734-2.0 Number of live births
% 2814-2.0 Ever used hormone-replacement therapy (HRT, 0 - no, 1 - yes, -1/-3 no answer))
%20116-2.0 Smoking status (-3 - no answer, 0 - never, 1 - previous, 2 - current)
% 1757-2.0 Facial ageing
%  991-2.0 Frequency of strenuous sports in last 4 weeks
% 1001-2.0 Duration of strenuous sports
% 1558-2.0 Alcohol intake frequency
% 1568-2.0 Average weekly red wine intake
% 1677-2.0 Breastfed as a baby
%   52-0.0 Month of birth
%   48-2.0 Waist circumference
%   49-2.0 Hip circumference
%   50-2.0 Standing height
% 4526-2.0 Happiness

codes_lifestyle = {'1070-2.0','1160-2.0','1239-2.0','1249-2.0','1289-2.0','1299-2.0','1309-2.0',...
                   '1319-2.0','1329-2.0','1349-2.0','1369-2.0','1379-2.0','1389-2.0','1558-2.0',...
                   '884-2.0','904-2.0'};
names_lifestyle = {'Time_spend_watching_TV','Sleep_duration','Current_tobacco_smoking','Past_tobacco_smoking',...
                   'Cooked_vegetable_intake','Salad_rawVegetable_intake','Fresh_fruit_intake',...
                   'Dried_fruit_intake','Oily_fish_intake','Processed_meat_intake','Beef_intake','Lamb_mutton_intake','Pork_intake',...
                   'Alcohol_intake_frequency','Moderate_activity','Vigorous_activity'};

codes_menopause = {'2814-2.0','2734-2.0','21001-2.0','3591-2.0','2834-2.0','6138-2.0','4079-2.0','4080-2.0','2443-2.0','738-2.0'};
names_menopause = {'HRT','Number of live births','BMI','Hysterectomy','Oophorectomy','Education','BP_diastolic','BP_systolic','Diabetes diagnosed by doctor','Income'};

switch abs(sel)
  case 0,  codes = {'3581-2.0','2714-2.0','2814-2.0',codes_menopause{:},codes_lifestyle{:}}; % analyzing variables in sample
           name = {'Analyzing sample',}; 
           name_codes = {'Age at Menopause','Age of Menarche [years]',names_menopause{:},'General Lifestyle'};
  case 1,  codes = {'3581-2.0','2714-2.0','2814-2.0'}; name = {'Age at Menopause [years]'}; % use all codes to obtain the same 1006 subjects for all analyses
  case 2,  codes = {'3581-2.0','2714-2.0','2814-2.0'}; name = {'Age of Menarche [years]'}; % use all codes to obtain the same 1006 subjects for all analyses
  case 3,  codes = {'3581-2.0','2714-2.0','2814-2.0'}; name = {'Reproductive Span [years]'}; % use all codes to obtain the same 1006 subjects for all analyses
  case 4,  codes = {'3581-2.0','2714-2.0','2814-2.0',codes_menopause{:},codes_lifestyle{:}}; 
           name = {'Age at Menopause [years]','Nuisance'}; % use all codes to obtain the same 1006 subjects for all analyses
  case 5,  codes = {'3581-2.0','2714-2.0','2814-2.0',codes_menopause{:},codes_lifestyle{:}}; 
           name = {'Age of Menarche [years]','Nuisance'}; % use all codes to obtain the same 1006 subjects for all analyses
  case 6,  codes = {'3581-2.0','2714-2.0','2814-2.0',codes_menopause{:},codes_lifestyle{:}}; 
           name = {'Reproductive Span [years]','Nuisance'}; % use all codes to obtain the same 1006 subjects for all analyses
  case 7,  codes = {codes_lifestyle{:}}; name = {'Lifestyle'};
  case 8,  codes = {'991-2.0'};  name = {'Frequency of strenuous sports in last 4 weeks'};
  case 9,  codes = {'1001-2.0'}; name = {'Duration of strenuous sports'};
  case 10, codes = {'991-2.0','1001-2.0'}; name = {'Hours per month of strenuous sports'};
  case 11, codes = {'4526-2.0'}; name = {'Happiness'};
  case 12, codes = {'48-2.0','49-2.0'}; name = {'Waist-hip ratio'};
  case 13, codes = {'52-0.0'}; name = {'Month of birth'};
  case 14, codes = {'1677-0.0'}; name = {'Breastfed as a baby'};
  case 15, codes = {'2814-2.0'}; name = {'Ever used hormone-replacement therapy'};
  case 16, codes = {'20116-0.0'}; name = {'Smoking status'};
  case 17, codes = {'2734-2.0'}; name = {'Number of live births'};
  case 18, codes = {'1757-2.0'}; name = {'Facial ageing'}; % 1 - younger, 2 - your age, 3 - older
end

% give warning that we already use that data for Eileen
if 0 && (abs(sel) < 8 || abs(sel) == 18)
  spm('alert','These are analyses that I prepared for Eileen and are only for demonstration')
end

num_in = [];
n = numel(codes);
ind_all = zeros(1,numel(header),'logical');

% go through all codes
for i = 1:n
  % find codes in header
  ind_codes = ismember(header,codes{i});
  % create array with numbers and get index where codes were found
  num_in = [num_in num(:,ind_codes)];
  ind_all = ind_all | ind_codes;
end

% negative entries means missing values, thus we set these to NaN
num_in(num_in<0) = NaN;

% do imputation for female data with menopause only
if abs(sel) < 7 && (abs(sel) > 3 || sel == 0)
  % get information about sex
  ind_codes = ismember(header,'31-0.0');
  female = ~num(:,ind_codes);

  fprintf('Impute data for %d females\n',sum(female))
  
  % ignore the first two columns with covariates of interest that should be not imputed
  num_in(female,4:end) = round(fillmissing(num_in(female,4:end),'knn'));
  fprintf('Code\t\tvalid data\timputations\n');
  for i=4:29
    fprintf('%10s\t%d\t%d\n',codes{i},sum(isfinite(num_in(female,i))),sum(female)-sum(isfinite(num_in(female,i))));
  end

  % single confounds that should be removed
  num_confound1 = num_in(:,4:numel(codes_menopause)+3);
  % combined lifestyle factor that should be removed
  num_confound2 = calculate_lifestyle_factor(num_in(:,(numel(codes_menopause)+4):end), names_lifestyle);
end

% get selected names from header
names = header(ind_all);

if 0 && numel(names) > n_covariates
  names = names(1:n_covariates);
end
  
% prepare some parameters for selected codes
switch abs(sel)
  case 0,  num_in1 = num_in(:,1); num_in1(num_in1<45 | num_in1 > 60) = NaN;
           num_in2 = num_in(:,2); num_in2(num_in2<10 | num_in2 > 18) = NaN;
           num_in3 = num_in(:,3); num_in3(num_in3<0) = NaN;
           num_in1(isnan(num_in2) | isnan(num_in3)) = NaN;
           num_in2(isnan(num_in1) | isnan(num_in3)) = NaN;
           num_in = [num_in1 - num_in2 num_in2  num_confound1 num_confound2];
  case 1,  num_in1 = num_in(:,1); num_in1(num_in1<45 | num_in1 > 60) = NaN;
           num_in2 = num_in(:,2); num_in2(num_in2<10 | num_in2 > 18) = NaN;
           num_in3 = num_in(:,3); num_in3(num_in3<0) = NaN;
           num_in1(isnan(num_in2) | isnan(num_in3)) = NaN;
           num_in = num_in1;
  case 2,  num_in1 = num_in(:,1); num_in1(num_in1<45 | num_in1 > 60) = NaN;
           num_in2 = num_in(:,2); num_in2(num_in2<10 | num_in2 > 18) = NaN;
           num_in3 = num_in(:,3); num_in3(num_in3<0) = NaN;
           num_in2(isnan(num_in1) | isnan(num_in3)) = NaN;
           num_in = num_in2;
  case 3,  num_in1 = num_in(:,1); num_in1(num_in1<45 | num_in1 > 60) = NaN;
           num_in2 = num_in(:,2); num_in2(num_in2<10 | num_in2 > 18) = NaN;
           num_in3 = num_in(:,3); num_in3(num_in3<0) = NaN;
           num_in1(isnan(num_in2) | isnan(num_in3)) = NaN;
           num_in2(isnan(num_in1) | isnan(num_in3)) = NaN;
           num_in = num_in1 - num_in2;
  case 4,  num_in1 = num_in(:,1); num_in1(num_in1<45 | num_in1 > 60) = NaN;
           num_in2 = num_in(:,2); num_in2(num_in2<10 | num_in2 > 18) = NaN;
           num_in3 = num_in(:,3); num_in3(num_in3<0) = NaN;
           num_in1(isnan(num_in2) | isnan(num_in3)) = NaN;
           num_in = [num_in1 num_confound1 num_confound2];
  case 5,  num_in1 = num_in(:,1); num_in1(num_in1<45 | num_in1 > 60) = NaN;
           num_in2 = num_in(:,2); num_in2(num_in2<10 | num_in2 > 18) = NaN;
           num_in3 = num_in(:,3); num_in3(num_in3<0) = NaN;
           num_in2(isnan(num_in1) | isnan(num_in3)) = NaN;
           num_in = [num_in2 num_confound1 num_confound2];
  case 6,  num_in1 = num_in(:,1); num_in1(num_in1<45 | num_in1 > 60) = NaN;
           num_in2 = num_in(:,2); num_in2(num_in2<10 | num_in2 > 18) = NaN;
           num_in3 = num_in(:,3); num_in3(num_in3<0) = NaN;
           num_in1(isnan(num_in2) | isnan(num_in3)) = NaN;
           num_in2(isnan(num_in1) | isnan(num_in3)) = NaN;
           num_in = [num_in1 - num_in2 num_confound1 num_confound2];
  case 7,  num_in = calculate_lifestyle_factor(num_in, names_lifestyle); order_poly = 1;
  case 8,  order_poly = 2;
  case 9,  num_in(num_in < 2) = NaN; order_poly = 2;
  case 10, num_in1 = num_in(:,1); num_in2 = num_in(:,2);
           freq = [1 2.5 4 10 18 30]; % frequency per month
           dur = [0.1 0.375 0.75 1.25 1.75 2.5 4]; % duration in hours
           % replace categories with frequency or duration
           for i=1:numel(freq)
             num_in1(num_in1 == i) = freq(i);
           end
           for i=1:numel(dur)
             num_in2(num_in2 == i) = dur(i);
           end
           % combine both by multiplying
           num_in = num_in1.*num_in2; 
  case 11, add_nuisance = 3;
  case 12, num_in1 = num_in(:,1);
           num_in2 = num_in(:,2);
           num_in = num_in1./num_in2; add_nuisance = 3;
  case 13, order_poly = 3; 
  case 14, force_regression = 0;
  case 15, force_regression = 0;
  case 16, force_regression = 0; num_in = num_in > 0;
  case 17, num_in(num_in > 5) = NaN;
  case 18, num_in1 = num_in; num_in1(num_in==3) = 2; num_in1(num_in==2) = 3; num_in = num_in1;
end

fprintf('%s\n',strjoin(name));
[num_array, age, sex, BA, BA_weighted] = cg_analyze_ukb_long(num_in, analysis, n_covariates, add_nuisance, ...
     order_poly, names, BA_weighting, force_regression);

if sel == 0
  for i=1:size(num_array,2)
    y = num_array(:,i);
    fprintf('\n%s:\n',name_codes{i});
    uniq_values = unique(y);
    n_uniq = numel(uniq_values);
    if n_uniq < 9
      h = hist(y,uniq_values);
      for j=1:n_uniq
        fprintf('Value %d n=%d\n',uniq_values(j),h(j));
      end
      fprintf('\n');
    end
    fprintf('Range=%g..%g\tMean=%g\tMedian=%g\tSD=%g\n',min(y),max(y),mean(y),median(y),std(y));
  end
end

set(gca,'FontSize',20)

% change Tick labels
switch abs(sel)
  case 18,   xlim([0.5 3.5]); set(gca,'XTick',1:3,'XTickLabel',{'younger','same','older'});
  case 11,   xlim([0.5 6.5]); set(gca,'XTick',1:6,'XTickLabel',{'Extremely happy','Very happy',...
                 'Moderately happy','Moderately unhappy','Very unhappy','Extremely unhappy'});
  otherwise, xlabel(name{1});
end

% change y-label
if sel > 0
  str = '_BL';
  ylabel('BrainAGE BL [years]');
  ylim([-15 15]);
else
  str = '_FU-BL';
  ylabel('Δ BrainAGE FU-BL [years]');
  ylim([-10 10]);
end

% save png and fig file
saveas(gcf,[strrep(strjoin(name),' ','_') str '_GPR.png']);
saveas(gcf,[strrep(strjoin(name),' ','_') str '.fig']);

if ~nargout
  clear num_array age male BA BA_weighted
end

function Lifestyle_score = calculate_lifestyle_factor(in, names_lifestyle)
% Calculate combined lifestyle factor based on Lancet paper and Claudia Barths R-script:
% https://www.thelancet.com/journals/lanpub/article/PIIS2468-2667(18)30200-7

if size(in,2) ~= 16
  error('Number of lifestyle parameters should be 16');
end

in(in < 0) = 0;

for i=1:16
  eval([names_lifestyle{i} '= in(:,i);']);
end

Smoking_binary = Current_tobacco_smoking > 0;
Alcohol_binary = Alcohol_intake_frequency == 1;
TV_watching_binary = Time_spend_watching_TV > 3;
Sleep_duration_binary = Sleep_duration > 0 & (Sleep_duration < 7 | Sleep_duration > 9);
Fruit_vegetable_intake = Dried_fruit_intake + Fresh_fruit_intake + Salad_rawVegetable_intake + Cooked_vegetable_intake;
Fruit_vegetable_intake_gramm = Fruit_vegetable_intake*80;
Fruit_vegetable_intake_binary = Fruit_vegetable_intake_gramm < 400;
Oily_fish_intake_binary = Oily_fish_intake < 2;
Red_meat_intake_binary = Beef_intake > 3| Pork_intake > 3 | Lamb_mutton_intake > 3;
Processed_meat_intake_binary = Processed_meat_intake > 2;
imaging_physical_activity = Moderate_activity < 5 | Vigorous_activity < 3;

Lifestyle_score = imaging_physical_activity + Processed_meat_intake_binary + Red_meat_intake_binary + ...
  Oily_fish_intake_binary + Fruit_vegetable_intake_binary + Sleep_duration_binary + TV_watching_binary + ...
  Alcohol_binary + Smoking_binary;
