function inds=get_subpop_indices(pop,m0)
%   GET_SUBPOP_INDICES sample a population to produce a subpopulation with
%   mean m0
sw=nanmean(diff(sort(pop))); %noise std
% sw=0;
w=sw*randn(size(pop)); %added noise
pop1=pop-m0+w;
[P,I0]=sort(abs(pop1));
s=(pop1(I0)>=0); %idicators of positive side
I1=[1;find(diff(s)~=0)+1];
inds=I0(I1);
nanmean(pop(inds));
end