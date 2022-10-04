#a note o how I made the pedigree file for lep_MAP3

cp mapping.txt PED4F_part.txt

sed -i "s/ /\t/g" PED4F_part.txt
sed -i "s/PX506\tPX553\t/PX553\tPX506\tPX553\tPX506\tPX506\tPX553\tPX506\tPX553\tA1_P\tA1_M\tA2_P\tA2_M\tB1_P\tB1_M\tB2_P\tB2_M\t/" PED4F_part.txt

cat PED4F_part.txt |sed -E "s/_[A-Z]+[0-9]*(\t)*/\t/g" - > PED4F.txt
sed -i "s/PX553\tPX506\tPX553\tPX506\tPX506\tPX553\tPX506\tPX553\t/A1\tA1\tA2\tA2\tB1\tB1\tB2\tB2\t" PED4F.txt

cat PED4F_part.txt >> PED4F.txt

cat PED4F_part.txt |sed -E "s/PX553|PX506/0/g" - | \
sed -E "s/A[12]_[PM]/PX553/g" - |sed -E "s/B[12]_[PM]/PX506/g" - | \
sed -E "s/A1_[A-Z]+[0-9]+\t/A1_P\t/g" - |sed -E "s/A2_[A-Z]+[0-9]+\t/A2_P\t/g" - | \
sed -E "s/B1_[A-Z]+[0-9]+\t/B1_P\t/g" - | sed -E "s/B2_[A-Z]+[0-9]+\t/B2_P\t/g" - | \
sed -E "s/B2_H11/B2_P/" - >> PED4F.txt


cat PED4F_part.txt |sed -E "s/PX553|PX506/0/g" - | \
sed -E "s/A[12]_[PM]/PX506/g" - |sed -E "s/B[12]_[PM]/PX553/g" - | \
sed -E "s/A1_[A-Z]+[0-9]+\t/A1_M\t/g" - |sed -E "s/A2_[A-Z]+[0-9]+\t/A2_M\t/g" - | \
sed -E "s/B1_[A-Z]+[0-9]+\t/B1_M\t/g" - | sed -E "s/B2_[A-Z]+[0-9]+\t/B2_M\t/g" - | \
sed -E "s/B2_H11/B2_M/" - >> PED4F.txt



cat PED4F_part.txt | \
sed -E "s/PX553\tPX506\tPX553\tPX506\tPX506\tPX553\tPX506\tPX553\tA1_P\tA1_M\tA2_P\tA2_M\tB1_P\tB1_M\tB2_P\tB2_M/1\t2\t1\t2\t1\t2\t1\t2\t1\t2\t1\t2\t1\t2\t1\t2/" - | \
sed -E "s/[AB12]+_[A-Z]+[0-9]+/2/g" - >> PED4F.txt


cat PED4F_part.txt | \
sed -E "s/PX553|PX506|[AB12]+_[PM]+/0/g" - | \
sed -E "s/[AB12]+_[A-Z]+[0-9]+/0/g" - >> PED4F.txt


sed -i "s/^/CHR\tPOS\t/g" PED4F.txt
