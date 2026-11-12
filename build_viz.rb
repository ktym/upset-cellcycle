#!/usr/bin/env ruby
# frozen_string_literal: true

# Build the dataset for the cell-cycle timeline viewer:
#   - reads the baseline TSV (61 curated interactions)
#   - adds the new candidate interactions from the survey
#   - attaches phase, human counterparts and disease notes
#   - writes candidates TSV (same columns as the baseline) and embeds JSON in the HTML template
require 'json'

OUT  = __dir__
BASE = "#{OUT}/sce_cellcycle_interactions.tsv"

PHASES = ['M/G1〜G1', 'Start（G1/S）', 'S期', 'G2/M', '中期（SAC）', '後期', '分裂終期'].freeze

# phase (start..end index into PHASES) for baseline rows, by numeric id
BASE_PHASE = {
  1..4 => [0, 0], 5..14 => [1, 1], 15..18 => [0, 0], 19..19 => [1, 1], 20..20 => [0, 0],
  21..21 => [1, 1], 22..22 => [0, 0], 23..26 => [2, 2], 27..33 => [3, 3], 34..36 => [4, 4],
  37..37 => [5, 5], 38..38 => [4, 4], 39..40 => [5, 5], 41..43 => [4, 4], 44..46 => [2, 2],
  47..51 => [6, 6], 52..52 => [4, 5], 53..56 => [6, 6], 57..60 => [5, 5], 61..61 => [0, 0]
}.freeze

# ORFs of genes that first appear in the candidates (from memory: verify against SGD)
NEW_ORF = {
  'WHI7' => 'YKR091W', 'NRM1' => 'YNR009W', 'STB1' => 'YNL309W', 'MSA1' => 'YOR066W', 'MSA2' => 'YKR077W',
  'YDJ1' => 'YNL064C', 'SSA1' => 'YAL005C', 'CDC48' => 'YDL126C', 'CKS1' => 'YBR135W', 'CSA1' => 'YPR174C',
  'CLB3' => 'YDL155W', 'MCM4' => 'YPR019W', 'MCM6' => 'YGL201C', 'SLD7' => 'YOR060C', 'RIF1' => 'YBR275C',
  'GLC7' => 'YER133W', 'MCM10' => 'YIL150C', 'MIH1' => 'YMR036C', 'ZDS1' => 'YMR273C', 'ZDS2' => 'YML109W',
  'RIM15' => 'YFL033C', 'IGO1' => 'YNL157W', 'IGO2' => 'YHR132W-A', 'MPS1' => 'YDL028C', 'SPC105' => 'YGL093W',
  'BUB1' => 'YGR188C', 'IPL1' => 'YPL209C', 'SLI15' => 'YBR156C', 'SGO1' => 'YOR073W', 'RTS1' => 'YOR014W',
  'KIN4' => 'YOR233W', 'LTE1' => 'YAL024C', 'BUD14' => 'YAR014C', 'NUD1' => 'YOR373W', 'AMN1' => 'YBR158W',
  'ACM1' => 'YPL267W', 'BMH1' => 'YER177W', 'BMH2' => 'YDR099W', 'HOG1' => 'YLR113W'
}.freeze

def ref(label, pmid: nil, url: nil) = { label: label, pmid: pmid, url: url }

# [module, sources, targets, type, 説明, refs, priority, [phase from, to]]
NEW = [
  ['G1/Start', %w[CLN3 CDC28], %w[WHI7], 'phosphorylation/inhibition',
   'Cln3-Cdc28がWhi7をリン酸化して不活性化する。Whi7はWhi5のパラログ。',
   [ref('Gomar-Alba et al. 2017', pmid: '28839131')], 'A', [0, 0]],
  ['G1/Start', %w[WHI7], %w[SWI4 SWI6], 'inhibition (binding)',
   'Whi7がSBFに結合してG1/S転写を抑える。Whi5と並列のStart抑制因子。',
   [ref('Gomar-Alba et al. 2017', pmid: '28839131')], 'A', [0, 0]],
  ['G1/Start', %w[PHO85], %w[WHI7], 'regulation (phosphorylation)',
   'Pho85によるWhi7の制御。作用の向きと部位は原著で要確認。',
   [ref('Ros-Carrero et al. 2024 EMBO Rep', url: 'https://doi.org/10.15252/embr.202357387')], 'B', [0, 0]],
  ['G1/Start', %w[YDJ1 SSA1], %w[CLN3], 'activation (release)',
   'シャペロンYdj1/Ssa1がCln3の放出と核移行を助ける。保持・放出の分子実体は見直しが進行中。',
   [ref('Vergés et al. 2007 Mol Cell 26:649'),
    ref('Aldea研 2019', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6467244/')], 'B', [0, 0]],
  ['G1/Start', %w[CDC48], %w[CLN3], 'activation (release)',
   'Cdc48がCln3の放出に関わる。',
   [ref('Parisi et al. 2018', pmid: '29950310')], 'B', [0, 0]],
  ['G1/Start (nutrient)', %w[CDC55], %w[WHI5], 'dephosphorylation/activation',
   'PP2A-Cdc55がWhi5を脱リン酸化してStartを遅らせる。貧栄養ではRim15-Igo1/2がこれを抑える。',
   [ref('Talarek et al. 2017 eLife 6:e26233', url: 'https://elifesciences.org/articles/26233')], 'B', [0, 0]],
  ['G1/Start (nutrient)', %w[CLN3 CDC28], %w[RIM15], 'inhibition',
   'Cln3-Cdk1がRim15を抑制し、Whi5リン酸化の恒常性を保つ。',
   [ref('Talarek et al. 2017 eLife 6:e26233', url: 'https://elifesciences.org/articles/26233')], 'B', [0, 0]],
  ['Stress (HOG)', %w[HOG1], %w[SIC1], 'phosphorylation (stabilization)',
   '浸透圧ストレスでHog1がSic1のT173をリン酸化し、ユビキチン化を妨げて安定化。G1で停止させる。',
   [ref('Escoté et al. 2004 Nat Cell Biol'), ref('Adrover et al. 2011 Sci Signal 4:ra63')], 'A', [0, 1]],
  ['G1/Start', %w[STB1], %w[SWI6], 'binding',
   'Stb1がSwi6を介してSBF/MBFに結合し、G1/S転写を調節する。',
   [ref('de Bruin et al. 2008', pmid: '18794370')], 'B', [1, 1]],
  ['G1/Start', %w[MSA1 MSA2], %w[SWI4 MBP1], 'binding/activation',
   'Msa1/Msa2がSBF・MBFに結合してG1/S転写のタイミングを調節する。',
   [ref('Ashe et al. 2008', pmid: '18160399'), ref('Miles et al. 2016', pmid: '27272642')], 'B', [1, 1]],
  ['CDK docking', %w[CKS1], %w[CDC28], 'binding',
   'Cks1がCdc28に結合し、リン酸化済み部位を足場に多点リン酸化を進める。Sic1分解のスイッチ性を担う。',
   [ref('Kõivomägi et al. 2011', pmid: '21993622')], 'A', [1, 1]],
  ['CDK docking', %w[CDC28 CLB5], %w[SIC1], 'phosphorylation',
   'Clb5-Cdc28がRxLドッキングを使ってSic1をリン酸化し、自らの阻害因子の分解を完成させる。',
   [ref('Kõivomägi et al. 2011', pmid: '21993622')], 'A', [1, 2]],
  ['S phase', %w[RIF1 GLC7], %w[MCM4], 'dephosphorylation/inhibition',
   'Rif1がPP1（Glc7）を呼び込み、DDKによるMcm4のリン酸化を打ち消して起点発火を抑える。',
   [ref('Hiraga et al. 2014 Genes Dev', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3937515/'),
    ref('Mattarocci et al. 2014', pmid: '24685139')], 'A', [1, 2]],
  ['S phase', %w[CDC7 DBF4], %w[MCM4 MCM6], 'phosphorylation',
   'DDKの必須標的はMcm4・Mcm6。再構成系で確定した。',
   [ref('Yeeles et al. 2015', pmid: '25739503')], 'A', [2, 2]],
  ['S phase', %w[SLD3 SLD7], %w[MCM4 MCM6], 'binding (recruitment)',
   'Sld3-Sld7がDDKでリン酸化されたMCMを読み取って結合し、Cdc45を呼び込む。',
   [ref('Deegan et al. 2016 EMBO J 35:961')], 'A', [2, 2]],
  ['S phase', %w[MCM10], %w[MCM2], 'binding/activation',
   'Mcm10がMCM複合体に直接結合し、CMGヘリカーゼの活性化に働く。',
   [ref('Douglas & Diffley 2016', pmid: '26719337')], 'B', [2, 2]],
  ['S phase checkpoint', %w[RAD53], %w[CDC7 DBF4], 'inhibition (binding block)',
   'Rad53がDDKのMCMへの結合を妨げる。',
   [ref('eLife 58571', url: 'https://elifesciences.org/articles/58571')], 'B', [2, 2]],
  ['S phase', %w[NRM1], %w[MBP1 SWI6], 'inhibition (binding)',
   'MBF標的であるNrm1がMBFに結合して転写を止める負のフィードバック。',
   [ref('de Bruin et al. 2006', pmid: '16916637')], 'A', [2, 2]],
  ['S phase checkpoint', %w[RAD53], %w[NRM1], 'phosphorylation/inhibition',
   '複製ストレス時にRad53がNrm1をリン酸化してMBFから外し、G1/S転写を維持する。',
   [ref('Travesa et al. 2012 EMBO J', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3321207/')], 'A', [2, 2]],
  ['APC/C', %w[ACM1 BMH1 BMH2], %w[CDH1], 'inhibition (binding)',
   'Acm1が14-3-3（Bmh1/2）と複合体を作り、KEN box・D boxを使う擬似基質としてCdh1を阻害する。G1後期からM期後半まで存在。',
   [ref('Martinez et al. 2006', pmid: '17030612'), ref('Enquist-Newman et al. 2008', pmid: '18498748'),
    ref('Burton et al. 2011 EMBO J 30:1818')], 'A', [1, 5]],
  ['APC/C', %w[CDC28], %w[ACM1], 'phosphorylation (stabilization)',
   'Cdc28によるリン酸化がBmh結合とAcm1の安定化に必要。',
   [ref('Ostapenko et al. 2008', pmid: '18519589')], 'A', [1, 4]],
  ['G2/M', %w[CSA1], %w[CLB3], 'binding (docking)',
   'Clb3特異的なPxxPxFドッキングモチーフをもつ基質Csa1（Ypr174c）。',
   [ref('Örd et al. 2020 Cell Rep', pmid: '32553169')], 'B', [3, 3]],
  ['G2/M', %w[MIH1], %w[CDC28], 'dephosphorylation/activation',
   'Mih1（Cdc25）がSwe1によるCdc28-Y19の抑制性リン酸化を外す。',
   [ref('Pal et al. 2008', pmid: '18316413')], 'A', [3, 3]],
  ['G2/M', %w[CDC55], %w[SWE1], 'dephosphorylation',
   'PP2A-Cdc55がCdk1によるSwe1のリン酸化に拮抗し、Swe1の活性化を制御する。',
   [ref('Harvey et al. 2011 MBoC 22:3595')], 'A', [3, 3]],
  ['G2/M', %w[ZDS1 ZDS2], %w[CDC55], 'binding/regulation',
   'Zds1/Zds2がPP2A-Cdc55に結合してその局在と働きを調節し、分裂期進入を助ける。',
   [ref('PMID 18762578', pmid: '18762578'), ref('Calabria et al. 2012', pmid: '22427694')], 'B', [3, 3]],
  ['G2/M', %w[RIM15], %w[IGO1 IGO2], 'phosphorylation/activation',
   'Rim15（Greatwall）がIgo1/2（エンドスルフィン）をリン酸化する。',
   [ref('Juanes et al. 2013', pmid: '23861665')], 'B', [3, 3]],
  ['G2/M', %w[IGO1 IGO2], %w[CDC55], 'inhibition (binding)',
   'リン酸化Igo1/2がPP2A-Cdc55を阻害する。',
   [ref('Juanes et al. 2013', pmid: '23861665')], 'B', [3, 3]],
  ['Stress (HOG)', %w[HOG1], %w[HSL1], 'phosphorylation/inhibition',
   'Hog1がHsl1のHsl7結合部位をリン酸化し、Hsl7がセプチン環から外れてSwe1が蓄積。G2で停止させる。',
   [ref('Clotet et al. 2006 EMBO J 25:2338', url: 'https://doi.org/10.1038/sj.emboj.7601095')], 'A', [3, 3]],
  ['M (SAC)', %w[MPS1], %w[SPC105], 'phosphorylation/activation',
   '未接着キネトコアでMps1がSpc105（KNL1）のMELTモチーフをリン酸化する。',
   [ref('London et al. 2012 Curr Biol', pmid: '22521787')], 'A', [4, 4]],
  ['M (SAC)', %w[SPC105], %w[BUB3 BUB1], 'binding (recruitment)',
   'Bub3がリン酸化MELTを直接読み取り、Bub3-Bub1がキネトコアに載る。',
   [ref('Primorac et al. 2013 eLife 2:e01030', url: 'https://elifesciences.org/articles/01030')], 'A', [4, 4]],
  ['M (SAC)', %w[MPS1], %w[BUB1], 'phosphorylation/activation',
   'キネトコア上のBub1もMps1でリン酸化され、Mad1-Mad2が結合できるようになる。',
   [ref('London & Biggins 2014; Moyle et al. 2014（Aravamudhan et al. 2016での引用）',
        url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5221577/')], 'B', [4, 4]],
  ['M (SAC)', %w[BUB1], %w[MAD1], 'binding (recruitment)',
   'リン酸化Bub1がMad1-Mad2を呼び込む。',
   [ref('London & Biggins 2014; Moyle et al. 2014（Aravamudhan et al. 2016での引用）',
        url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5221577/')], 'B', [4, 4]],
  ['M (SAC)', %w[IPL1 SLI15], %w[MAD2], 'activation (indirect)',
   'Ipl1/Aurora Bが張力のない接着を不安定化し、SACを点灯させる。キネトコア経由の間接効果。',
   [ref('Biggins & Murray 2001 Genes Dev 15:3118')], 'B', [4, 4]],
  ['M (SAC)', %w[BUB1], %w[SGO1 RTS1], 'activation (recruitment)',
   'Bub1がH2A-S121をリン酸化してSgo1-PP2A(Rts1)をペリセントロメアに呼ぶ。張力が生じると外れる。',
   [ref('Nerusheva et al. 2014 Genes Dev 28:1291', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4066400/'),
    ref('Genetics 2017 205:1169', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5340331/')], 'B', [4, 4]],
  ['M (SAC)', %w[GLC7], %w[SPC105], 'dephosphorylation/inhibition',
   'PP1（Glc7）がMELTのリン酸化を外し、SACを解除する。',
   [ref('London et al. 2012 Curr Biol', pmid: '22521787')], 'A', [5, 5]],
  ['APC/C', %w[CDC20], %w[ACM1], 'ubiquitination/proteolysis',
   'APC/C-Cdc20が後期にAcm1を分解し、Cdh1の活性化を順序づける。',
   [ref('Enquist-Newman et al. 2008', pmid: '18498748')], 'A', [5, 5]],
  ['SPOC', %w[KIN4], %w[BFA1], 'phosphorylation/activation',
   '紡錘体位置がずれるとKin4がBfa1をリン酸化し、Cdc5による抑制性リン酸化を受けないようにする。',
   [ref('Maekawa et al. 2007; Caydasi & Pereira 2009（eLife 72833の整理）',
        url: 'https://elifesciences.org/articles/72833')], 'B', [5, 5]],
  ['SPOC', %w[LTE1], %w[KIN4], 'inhibition',
   '芽側でLte1がKin4を娘側SPBから遠ざけ、Tem1を活性化させる。',
   [ref('Bertazzi et al. 2011; Falk et al. 2011（総説での引用）',
        url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7185974/')], 'B', [5, 5]],
  ['SPOC', %w[GLC7 BUD14], %w[BFA1], 'dephosphorylation/activation',
   'Glc7-Bud14がBfa1を脱リン酸化してSPOCのエフェクターを完全に活性化する。',
   [ref('eLife 72833', url: 'https://elifesciences.org/articles/72833')], 'B', [5, 5]],
  ['Mitotic exit (MEN)', %w[CDC15], %w[NUD1], 'phosphorylation',
   'Cdc15がSPBスキャフォールドNud1をリン酸化し、Dbf2-Mob1のドッキング部位を作る。',
   [ref('Rock et al. 2013（総説での引用）', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7185974/')], 'A', [6, 6]],
  ['Mitotic exit (MEN)', %w[NUD1], %w[DBF2 MOB1], 'binding (recruitment)',
   'リン酸化Nud1がDbf2-Mob1をSPBに呼び、Cdc15による活性化を可能にする。',
   [ref('Rock et al. 2013（総説での引用）', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7185974/')], 'A', [6, 6]],
  ['Mitotic exit (MEN)', %w[DBF2 MOB1], %w[CDC14], 'phosphorylation',
   'Dbf2-Mob1がCdc14のNLS近傍をリン酸化してNLSを無効化し、細胞質に留める。',
   [ref('Mohl et al. 2009 JCB 184:527', pmid: '19221193')], 'A', [6, 6]],
  ['Mitotic exit (MEN)', %w[DBF2 MOB1], %w[NET1], 'phosphorylation',
   'Cdc5のプライミング後、核小体に入ったDbf2-Mob1がNet1をリン酸化してCdc14を放出させる。',
   [ref('Zhou et al. 2021', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7822594/')], 'B', [6, 6]],
  ['Mitotic exit (reset)', %w[AMN1], %w[TEM1], 'inhibition (binding, degradation)',
   '娘細胞特異的に発現するAmn1がTem1に結合してCdc15との結合を妨げ、SPBから外して分解にも導く。',
   [ref('Wang et al. 2003 Cell 112:697', pmid: '12628189'),
    ref('Devault & Piatti 2021 J Cell Sci', pmid: '34518877')], 'B', [6, 6]]
].freeze

# Human counterparts (from memory: to be replaced by Alliance / Monarch orthology calls)
#   o: orthologue, f: functional counterpart without clear sequence orthology
HUMAN = {
  'CDC28' => ['o', %w[CDK1 CDK2]], 'CLN1' => ['f', %w[CCNE1]], 'CLN2' => ['f', %w[CCNE1]], 'CLN3' => ['f', %w[CCND1]],
  'CLB5' => ['f', %w[CCNA2 CCNE1]], 'CLB6' => ['f', %w[CCNA2]], 'CLB2' => ['o', %w[CCNB1]], 'CLB3' => ['o', %w[CCNB1]],
  'WHI5' => ['f', %w[RB1]], 'WHI7' => ['f', %w[RB1]], 'SWI4' => ['f', %w[E2F1]], 'SWI6' => ['f', %w[E2F1]],
  'MBP1' => ['f', %w[E2F1]], 'SIC1' => ['f', %w[CDKN1B]], 'CKS1' => ['o', %w[CKS1B]], 'CAK1' => ['f', %w[CDK7]],
  'CDC4' => ['o', %w[FBXW7]], 'SKP1' => ['o', %w[SKP1]], 'CDC53' => ['o', %w[CUL1]], 'CDC34' => ['o', %w[CDC34]],
  'HRT1' => ['o', %w[RBX1]], 'PHO85' => ['o', %w[CDK5]], 'CDC6' => ['o', %w[CDC6]], 'ORC1' => ['o', %w[ORC1]],
  'TAH11' => ['o', %w[CDT1]], 'MCM2' => ['o', %w[MCM2]], 'MCM4' => ['o', %w[MCM4]], 'MCM6' => ['o', %w[MCM6]],
  'MCM10' => ['o', %w[MCM10]], 'CDC7' => ['o', %w[CDC7]], 'DBF4' => ['o', %w[DBF4]], 'SLD2' => ['o', %w[RECQL4]],
  'SLD3' => ['o', %w[TICRR]], 'SLD7' => ['o', %w[MTBP]], 'RIF1' => ['o', %w[RIF1]], 'GLC7' => ['o', %w[PPP1CA]],
  'RAD53' => ['o', %w[CHEK2]], 'CHK1' => ['o', %w[CHEK1]], 'MEC1' => ['o', %w[ATR]], 'TEL1' => ['o', %w[ATM]],
  'RAD9' => ['f', %w[TP53BP1]], 'MRC1' => ['o', %w[CLSPN]], 'MCM1' => ['o', %w[SRF]], 'FKH2' => ['f', %w[FOXM1]],
  'SWE1' => ['o', %w[WEE1]], 'MIH1' => ['o', %w[CDC25C]], 'HSL7' => ['o', %w[PRMT5]], 'CDC5' => ['o', %w[PLK1]],
  'CDC55' => ['o', %w[PPP2R2A]], 'RTS1' => ['o', %w[PPP2R5D]], 'RIM15' => ['o', %w[MASTL]],
  'IGO1' => ['o', %w[ENSA]], 'IGO2' => ['o', %w[ENSA]], 'MAD1' => ['o', %w[MAD1L1]], 'MAD2' => ['o', %w[MAD2L1]],
  'MAD3' => ['o', %w[BUB1B]], 'BUB1' => ['o', %w[BUB1]], 'BUB3' => ['o', %w[BUB3]], 'MPS1' => ['o', %w[TTK]],
  'SPC105' => ['o', %w[KNL1]], 'IPL1' => ['o', %w[AURKB]], 'SLI15' => ['o', %w[INCENP]], 'SGO1' => ['o', %w[SGO1]],
  'CDC20' => ['o', %w[CDC20]], 'CDH1' => ['o', %w[FZR1]], 'CDC16' => ['o', %w[CDC16]], 'CDC23' => ['o', %w[CDC23]],
  'CDC27' => ['o', %w[CDC27]], 'PDS1' => ['f', %w[PTTG1]], 'ESP1' => ['o', %w[ESPL1]], 'MCD1' => ['o', %w[RAD21]],
  'CDC14' => ['o', %w[CDC14A]], 'CDC15' => ['f', %w[STK4]], 'DBF2' => ['o', %w[LATS1]], 'MOB1' => ['o', %w[MOB1A]],
  'NUD1' => ['f', %w[CNTRL]], 'ACM1' => ['f', %w[FBXO5]], 'BMH1' => ['o', %w[YWHAE]], 'BMH2' => ['o', %w[YWHAE]],
  'HOG1' => ['o', %w[MAPK14]], 'YDJ1' => ['o', %w[DNAJA1]], 'SSA1' => ['o', %w[HSPA1A]], 'CDC48' => ['o', %w[VCP]]
}.freeze

# Disease / clinical notes for human genes (from memory: replace with Monarch / ClinVar / TogoVar lookups)
DISEASE = {
  'RB1' => '網膜芽細胞腫。多くのがんで不活性化', 'CCND1' => '多くのがんで増幅・転座（マントル細胞リンパ腫など）',
  'CCNE1' => '卵巣がんなどで増幅', 'CDKN1B' => '多発性内分泌腫瘍症4型（MEN4）',
  'FBXW7' => 'がん抑制遺伝子（T-ALL、大腸がんなどで変異）', 'CKS1B' => '多発性骨髄腫の1q増幅領域',
  'BUB1B' => 'mosaic variegated aneuploidy症候群1', 'KNL1' => '原発性小頭症（MCPH4）',
  'ATM' => '毛細血管拡張性運動失調症。乳がんなどの易罹患性', 'ATR' => 'Seckel症候群1',
  'CHEK2' => '乳がんなどの遺伝性易罹患性', 'RECQL4' => 'Rothmund-Thomson症候群',
  'ORC1' => 'Meier-Gorlin症候群', 'CDC6' => 'Meier-Gorlin症候群', 'CDT1' => 'Meier-Gorlin症候群',
  'MCM4' => 'NK細胞欠損と副腎不全を伴う免疫不全', 'RAD21' => 'Cornelia de Lange症候群4',
  'VCP' => '封入体ミオパチー・前頭側頭型認知症（IBMPFD）、ALS', 'CDC14A' => '常染色体劣性難聴（DFNB32）',
  'LATS1' => 'Hippo経路のがん抑制因子', 'STK4' => 'STK4欠損症（原発性免疫不全）',
  'PPP2R5D' => '知的障害を伴う神経発達症', 'PLK1' => '抗がん剤の標的（PLK1阻害薬）',
  'WEE1' => '抗がん剤の標的（WEE1阻害薬）', 'AURKB' => '抗がん剤の標的（Aurora阻害薬）',
  'CDC7' => '抗がん剤の標的（CDC7阻害薬）', 'TTK' => '抗がん剤の標的（MPS1阻害薬）',
  'CDK1' => 'CDK阻害薬の標的', 'CDK2' => 'CDK阻害薬の標的', 'CHEK1' => '抗がん剤の標的（CHK1阻害薬）',
  'PTTG1' => '下垂体腫瘍などで過剰発現'
}.freeze

def classify(type)
  case type
  when /inhibition|repression|ubiquitination|proteolysis/ then 'inh'
  when /activation|expression|release|recruitment|stabilization/ then 'act'
  when /phosph|regulation/ then 'mod'
  else 'bind'
  end
end

lines = File.readlines(BASE, chomp: true, encoding: 'UTF-8')
fields = lines.shift.split("\t", -1)
base = lines.map { |l| fields.zip(l.split("\t", -1)).to_h }

orf = {}
base.each do |r|
  %w[source target].each do |k|
    r["#{k}_genes"].split('|').zip(r["#{k}_orfs"].split('|')).each { |g, o| orf[g] = o }
  end
end
orf.merge!(NEW_ORF) { |_, a, _| a }

ix = base.map do |r|
  n = r['id'][/\d+\z/].to_i
  pmids = r['pmids'].split(';').reject(&:empty?)
  labels = r['references'].split(/;\s*/)
  refs = pmids.each_with_index.map { |p, i| ref(labels[i] || "PMID #{p}", pmid: p) }
  refs = labels.map { |l| ref(l) } if refs.empty?
  { id: r['id'], set: 'base', pri: nil, module: r['module'], src: r['source_genes'].split('|'),
    tgt: r['target_genes'].split('|'), type: r['interaction_type'], cls: classify(r['interaction_type']),
    desc: r['description_ja'], refs: refs, ph: BASE_PHASE.find { |k, _| k.cover?(n) }[1] }
end

cand_rows = []
NEW.each_with_index do |(mod, src, tgt, type, desc, refs, pri, ph), i|
  id = format('SCE-CC-N%02d', i + 1)
  (src + tgt).each { |g| raise "no ORF for #{g}" unless orf[g] }
  ix << { id: id, set: 'new', pri: pri, module: mod, src: src, tgt: tgt, type: type, cls: classify(type),
          desc: desc, refs: refs, ph: ph }
  cand_rows << [id, mod, src.join('|'), src.map { |g| orf[g] }.join('|'), tgt.join('|'),
                tgt.map { |g| orf[g] }.join('|'), type, "#{desc}（優先度#{pri}）",
                refs.filter_map { |r| r[:pmid] }.join(';'),
                refs.map { |r| r[:url] ? "#{r[:label]} <#{r[:url]}>" : r[:label] }.join('; '),
                'WebSearch（サーベイ時に確認）', '未照合', '未照合']
end

File.open("#{OUT}/sce_cellcycle_candidates.tsv", 'w', encoding: 'UTF-8') do |f|
  f.puts fields.join("\t")
  cand_rows.each { |r| f.puts r.join("\t") }
end

genes = ix.flat_map { |x| x[:src] + x[:tgt] }.uniq.to_h do |g|
  rel, hs = HUMAN[g]
  [g, { orf: orf[g], rel: rel, human: hs || [] }]
end
data = { phases: PHASES, genes: genes, disease: DISEASE, ix: ix }
html = File.read("#{OUT}/viewer_template.html", encoding: 'UTF-8')
File.write("#{OUT}/cellcycle_timeline.html", html.sub('/*DATA*/null') { JSON.generate(data) }, encoding: 'UTF-8')
warn "#{ix.size} interactions (#{cand_rows.size} new), #{genes.size} genes, #{genes.count { |_, v| v[:human].any? }} with human counterpart"
