function BA_menopause

xls_file = '/Users/gaser/Dropbox/UKBiobank/Jena_long3046/tables/ukb49261_long3046.xlsx';
[num, header] = cg_get_ukb_data(xls_file);

for i=[3 1 2 6 4 5]
  [num_array, age, sex, BA, BA_weighted] = cg_analyze_ukb_long_BAonly(i, num,header);
  [num_array, age, sex, BA, BA_weighted] = cg_analyze_ukb_long_BAonly(-i,num,header);
  cat_plot_scatter(num_array(:,1),BA_weighted,'Color',age,'plottype','scatter','jitter',1)
  colorbar
end
