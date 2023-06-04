function  [pk_ts,pk_amp,numpk,pval,pcnt,zext_ts,pk_ts2,pk_amp2,numpk2] = detect_RSSevents(ts,R,offsets,pthr,zext)

N = size(ts,2);
lts = size(ts,1);

% compute ets and rss
ets = fcn_edgets(ts); ets(isnan(ets))=0;
rssts = sum(ets.^2,2).^0.5;

% circshift null
pcnt = zeros(lts,1);
for r=1:R
    tsr = zeros(lts,N);
    for n=1:N
        tsr(:,n) = circshift(ts(:,n),offsets(randi(length(offsets))));
    end
    etsr = fcn_edgets(tsr); etsr(isnan(etsr))=0;
    rsstsr = sum(etsr.^2,2).^0.5;
    pcnt = pcnt+(rssts>max(rsstsr));        % '1' if rssts>all rsstsr values
end

% pval
pval = 1-pcnt./R;

% determine peaks as intersection of 'findpeaks' and 'pvals'
% findpeaks
[~, fp_ts] = findpeaks(rssts);
% intersection
pk_ts = intersect(find(pval<pthr),fp_ts);
pk_amp = rssts(pk_ts);
numpk = length(pk_ts);

% flag those 'events' that may be due to extreme z-scores (> abs(zthr))
zext_ts = union(find(sum(ts'>zext)),find(sum(ts'<-zext)));
mask1 = zeros(lts,1); mask2 = zeros(lts,1);
mask1(pk_ts) = 1; mask2(zext_ts) = 1;
pk_ts2 = find((mask1==1));%&(mask2==0));    % pick events that are not due to extreme z-scores
pk_amp2 = rssts(pk_ts2);
numpk2 = length(pk_ts2);

