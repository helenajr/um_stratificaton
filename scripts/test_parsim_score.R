# check that temp is always == 0 ---------------------------------------------------

# young and disomy ------------------------------------------------------------------
temp <- filter(data_lp, ageatpm < 81, chr3 == "Disomy", lbd < 10.1, score != 0.7)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Disomy", lbd >= 10.1 & lbd <= 12, score != 1.7)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Disomy", lbd >= 12.1 & lbd <= 14, score != 3.3)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Disomy", lbd >= 14.1 & lbd <= 16, score != 5.3)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Disomy", lbd >= 16.1 & lbd <= 18, score != 12)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Disomy", lbd >= 18.1 & lbd <= 28, score != 12)

# young and monosomy ------------------------------------------------------------------------------
temp <- filter(data_lp, ageatpm < 81, chr3 == "Monosomy", lbd < 10.1, score != 13)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Monosomy", lbd >= 10.1 & lbd <= 12, score != 20)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Monosomy", lbd >= 12.1 & lbd <= 14, score != 33)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Monosomy", lbd >= 14.1 & lbd <= 16, score != 42)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Monosomy", lbd >= 16.1 & lbd <= 18, score != 59)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Monosomy", lbd >= 18.1 & lbd <= 28, score != 70)

# young and unknown -------------------------------------------------------------------------------
temp <- filter(data_lp, ageatpm < 81, chr3 == "Unknown", lbd < 10.1, score != 3.1)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Unknown", lbd >= 10.1 & lbd <= 12, score != 7)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Unknown", lbd >= 12.1 & lbd <= 14, score != 15)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Unknown", lbd >= 14.1 & lbd <= 16, score != 25)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Unknown", lbd >= 16.1 & lbd <= 18, score != 43)
temp <- filter(data_lp, ageatpm < 81, chr3 == "Unknown", lbd >= 18.1 & lbd <= 28, score != 56)

# old and disomy ----------------------------------------------------------------------------------
temp <- filter(data_lp, ageatpm > 80, chr3 == "Disomy", lbd < 10.1, score != 0.6)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Disomy", lbd >= 10.1 & lbd <= 12, score != 1.5)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Disomy", lbd >= 12.1 & lbd <= 14, score != 2.9)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Disomy", lbd >= 14.1 & lbd <= 16, score != 4.6)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Disomy", lbd >= 16.1 & lbd <= 18, score != 11)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Disomy", lbd >= 18.1 & lbd <= 28, score != 11)

# old and monosomy
temp <- filter(data_lp, ageatpm > 80, chr3 == "Monosomy", lbd < 10.1, score != 11)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Monosomy", lbd >= 10.1 & lbd <= 12, score != 18)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Monosomy", lbd >= 12.1 & lbd <= 14, score != 29)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Monosomy", lbd >= 14.1 & lbd <= 16, score != 37)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Monosomy", lbd >= 16.1 & lbd <= 18, score != 53)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Monosomy", lbd >= 18.1 & lbd <= 28, score != 64)

# old and unknown
temp <- filter(data_lp, ageatpm > 80, chr3 == "Unknown", lbd < 10.1, score != 2.7)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Unknown", lbd >= 10.1 & lbd <= 12, score != 6.1)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Unknown", lbd >= 12.1 & lbd <= 14, score != 14)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Unknown", lbd >= 14.1 & lbd <= 16, score != 22)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Unknown", lbd >= 16.1 & lbd <= 18, score != 38)
temp <- filter(data_lp, ageatpm > 80, chr3 == "Unknown", lbd >= 18.1 & lbd <= 28, score != 51)

temp <- filter(data_lp, is.na(score))


