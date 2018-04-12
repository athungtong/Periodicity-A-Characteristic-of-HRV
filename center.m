function s=center(s)
% % s=center(s) normalized time series input s so that it has zero mean and
% unit variance.
s=(s-mean(s))/std(s);
