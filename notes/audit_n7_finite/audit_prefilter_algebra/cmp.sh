D=/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_prefilter_algebra
cd $D
cat n7U_indep_0*.txt | awk '{print $1, $2, $5}' > cur.txt
echo "done_lines $(wc -l < cur.txt)  pos $(awk '$2=="POS"' cur.txt | wc -l)  never $(awk '$2=="NEVER"' cur.txt | wc -l)  not_df1_drop0 $(cat n7U_indep_0*.txt | grep -v -c "df=1 drop=0$")"
# disagreements: mine POS but resolver dropped, or mine NEVER but resolver R
awk 'NR==FNR{r[$1]=1; next} ($2=="POS" && !($1 in r)) || ($2=="NEVER" && ($1 in r))' n7_R_resolver.txt cur.txt > disagree.txt
echo "disagreements $(wc -l < disagree.txt)"
