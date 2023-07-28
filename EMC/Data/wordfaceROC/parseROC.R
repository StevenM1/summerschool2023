# The S column is named as ‘type’, A named as ‘cond’, R named as ‘confresponse’ 
# in our files. The ratings goes from 1 to 6 (from definitely new to definitely old). 
# I also add the RT column there in case you need it.
# I excluded trials based on these 3 criteria:
# 1. Subjects with bad performance (too low d’ or too fast RTs)
# 2. RTs<0.3s and RTs>4s
# 3. ’99' trials
# 
# In type, 1 is a new stimulus and 2 is a old stimulus. I got the order of the 
# confresponse column wrong in previous email, we actually recorded them as 
# definitely old upwards (1 = def old and 6 is def new).


WordFace <- read.csv("Data/wordfaceROC.csv")
# length(unique(WordFace$subj)) # 54
# # remove non-responses
NR <- WordFace$confresponse==99
# round(sort(100*tapply(NR,WordFace$subj,mean)),1)
# # 104 106 107 109 116 117 120 121 126 128 130 131 134 136 143 144 147 148 153 156 157 161 101 102 122 
# # 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.1 0.1 
# # 129 132 135 138 112 115 125 142 155 137 141 110 149 113 139 146 154 151 123 119 124 
# # 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.3 0.3 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.9 2.3 2.7 
WordFace <- WordFace[!NR,]

# reverse confresponse
WordFace$confresponse <-  7-WordFace$confresponse
WordFace$response <-  WordFace$confresponse>3 # OLD RESPONSE
WordFace$C <- (WordFace$type==1 & !WordFace$response) | (WordFace$type==2 & WordFace$response)

# badHaomin <- c(103, 111, 118, 145, 127, 105, 140, 150)
# # Keep in single session, 103, 111, 118, 127, 145 were excluded for having 
# # <0.5 d’ in words and <0.2 d’ in faces, 105, 140, 150 didnt use confidence scale
# BUT LEAVES 9 WITH d'(faces) < .2

# Keep only sig better than chance on words and faces (always fails on faces)
nonsig <- tapply(WordFace$C,WordFace[,c("subj","cond")],function(x){binom.test(sum(x),length(x))$p.value > .05})
bad <- row.names(nonsig)[apply(nonsig,1,any)]
#  [1] "101" "102" "103" "110" "111" "113" "118" "123" "127" "132" "135" "139" "140"
# [14] "145" "150"
WordFace <- WordFace[!(WordFace$subj %in% bad),]

# # Remove non-compliant use of confidence scale
# table(WordFace$confresp,WordFace$subj,WordFace$cond)
# # , ,  = faces
# # 
# #    
# #     104 105 106 107 109 112 115 116 117 119 120 121 122 124 125 126 128 129 130 131
# #   1  20 351 244 412 390 213 263 289 152   0  66  48 301 158  41  39 167 179 121  52
# #   2  58   0 100  62   4   8  34  24 218   0 316 183  42 134  50 117  80 131 167 244
# #   3 146   0  49  17  22  93  27  19  47 192  77  49  72  91  65 391  55  58 114 148
# #   4  90   0  71  41  11  10   7  13  45 534 122  53  89  85  63 124 104  88 156 152
# #   5  42   0  72  23  16   6  48  15 103  19 137 200  82 129  84  64  82  78 124 115
# #   6  28 415 232 213 325 436 387 408 203   1  50 235 181 156  80  33 280 234  86  57
# #    
# #     134 136 137 138 141 142 143 144 146 147 148 149 151 153 154 155 156 157 161
# #   1 234 135 359   6  46 259  70  22 132 160  11  92  97 344 297 104 457  17  61
# #   2 169 181  47 113  99  16 143 151 126  86  52 108  46  69   6 227   9 234  77
# #   3  91  49  47 312 266  48 168   2  84  69 297  23  37   7  21  61   0  97  45
# #   4  42  69  21 179 237  93 185   9  92 117 161   0  54  27 239  74   0 113  39
# #   5 101 164  32  93  73  65  92  69 186 102 150  57 319  78   6 182  22 249  50
# #   6 131 170 257  65  45 286 110 131 145 234  97 483 209 243 198 120 280  58 112
# # 
# # , ,  = words
# # 
# #    
# #     104 105 106 107 109 112 115 116 117 119 120 121 122 124 125 126 128 129 130 131
# #   1  11 435 405 530 493 204 459 357 203   0 163 123 423 204  74  27 233 261 157  17
# #   2  29   0  28   9   0  10  27  18 201   0 262 156  11 134  64 149  65  82 147 289
# #   3 222   1  13   0  26  98  78   7  32 453  31  23  14  83  74 399  47  44 106 110
# #   4  61   1   9  22  12   6  28   5  19 282  95  53  19 100  40  71  46  39 138 142
# #   5  30   0  20   9   7  10  34   5 110  18  99  81  20  85  62  58  57  31 110 106
# #   6  31 331 293 198 230 440 142 376 203   1 118 332 281 136  70  64 320 310 110 104
# #    
# #     134 136 137 138 141 142 143 144 146 147 148 149 151 153 154 155 156 157 161
# #   1 329 158 481  15 188 253 220  15 118 209 111 341 145 445 303 134 495 312 128
# #   2 110 197  15 196  85  18 172 143 157  80  61  92 137   3   9 190   7  91  48
# #   3  30  17  14 313 157  63 138  11  94  51 263   7  72   0   3  48   1   6   7
# #   4  25  33  15 100  99  82 109  15  62  43  49   3  32   1 256  82   0   7  15
# #   5  68 119   3  24  54  45  45  85 157  23  72  41 164   7   9 100   5  55  27
# #   6 206 244 240 119 182 306  84 115 175 362 212 283 215 312 181 212 260 297 159

  
badConf <- c(105,119,156)
WordFace <- WordFace[!(WordFace$subj %in% badConf),]


# sort(table(WordFace$subj))  # 36 left
# #  125  104  144  161  124  151  146  154  149  137  141  112  115  142  155  122  129  138 
# #  767  768  768  768 1495 1527 1528 1528 1530 1531 1531 1534 1534 1534 1534 1535 1535 1535 
# #  106  107  109  116  117  120  121  126  128  130  131  134  136  143  147  148  153  157 
# # 1536 1536 1536 1536 1536 1536 1536 1536 1536 1536 1536 1536 1536 1536 1536 1536 1536 1536 


fast <- WordFace$RT<.3
# slow <- WordFace$RT>4

# hist(WordFace$RT,breaks="fd")
# par(mfrow=c(3,4))
# for (i in unique(WordFace$subj)) {
#    hist(WordFace$RT[WordFace$subj==i],breaks=seq(0,8,by=.1),main=i)
#   abline(v=.3)
#   abline(v=4)
# }

uc <- tapply(WordFace$RT,WordFace$subj,mean) + 3*tapply(WordFace$RT,WordFace$subj,sd)
slow <- factor(WordFace$subj)
levels(slow) <- uc[levels(slow)]
slow <- WordFace$RT > as.numeric(as.character(slow))


round(100*sort(tapply(fast,WordFace$subj,mean)),2)
#  104  107  109  115  116  117  120  121  122  126  128  129  130  131  134  136  138  141 
# 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 
#  142  143  144  147  148  151  153  155  157  161  106  112  149  154  125  137  146  124 
# 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.07 0.13 0.13 0.20 0.26 0.26 0.46 0.47  

round(100*sort(tapply(slow,WordFace$subj,mean)),2)
#  126  134  143  161  138  155  116  136  146  154  106  109  149  120  147  142  107  130  148  157  104 
# 1.11 1.50 1.50 1.56 1.56 1.63 1.69 1.69 1.77 1.77 1.82 1.82 1.83 1.89 1.89 1.96 2.02 2.02 2.02 2.02 2.08 
#  121  131  128  144  125  117  122  129  112  115  153  137  151  124  141 
# 2.08 2.08 2.21 2.21 2.22 2.28 2.28 2.28 2.35 2.35 2.41 2.48 2.49 2.68 3.53  

WordFace <- WordFace[!fast & !slow,]

dp <- apply(qnorm(tapply(WordFace$response,WordFace[,c("subj","cond","type")],mean)),1:2,diff)
round(sort(dp[,"words"]),2)
#  124  130  146  112  154  125  115  126  117  144  155  151  131  121  136  149  128  141  122  104  107 
# 0.56 0.57 0.70 0.75 0.81 0.93 0.95 0.99 1.01 1.05 1.12 1.15 1.26 1.27 1.31 1.33 1.45 1.46 1.52 1.53 1.67 
#  106  120  142  138  143  147  109  137  116  134  129  157  153  148  161 
# 1.68 1.70 1.72 1.73 1.80 1.89 1.96 2.05 2.08 2.33 2.38 2.65 2.70 2.81 2.85 
round(sort(dp[,"faces"]),2)
#  125  147  144  149  130  151  155  124  116  141  106  112  154  137  122  146  117  126  131  120  121 
# 0.27 0.28 0.29 0.31 0.31 0.32 0.34 0.35 0.35 0.35 0.35 0.36 0.39 0.40 0.41 0.41 0.42 0.49 0.50 0.51 0.55 
#  157  115  153  142  107  128  104  138  136  109  143  161  148  129  134 
# 0.58 0.73 0.76 0.78 0.82 0.83 0.91 0.92 0.94 0.98 0.98 0.99 1.07 1.19 1.40 

# Make into EMC format
wordfaceROC <- WordFace[,c(1,3:5,2)]
names(wordfaceROC) <- c("subjects","FW","S","R","rt")
wordfaceROC$subjects <- factor(wordfaceROC$subjects)
wordfaceROC$R <- factor(wordfaceROC$R)
wordfaceROC$S <- factor(wordfaceROC$S,labels=c("new","old"))
save(wordfaceROC,file="Data/wordfaceROC.RData")
