#!/usr/bin/env ruby
# frozen_string_literal: true

# Check the interactions of one master table against KEGG and Reactome and record the
# outcome, with the date, in the table's `evidence` column. The `curation` column
# (TODO / DONE, manual literature check) is never touched.
#
#   ruby check.rb cellcycle_hsa.tsv --kgml hsa04110.xml \
#        --reactome reactome.homo_sapiens.interactions.tab-delimited.txt --gene-list hsa.list
#   ruby check.rb cellcycle_sce.tsv --kgml sce04111.xml --kgml sce04110.xml \
#        --reactome reactome.all_species.interactions.tab-delimited.txt
#
# Inputs
#   TSV          cellcycle_<org>.tsv; <org> (KEGG organism code) is taken from the file name.
#   --kgml       KGML file of a pathway of that organism (repeatable):
#                  https://rest.kegg.jp/get/hsa04110/kgml
#   --reactome   Reactome interaction file, tab-delimited flavour:
#                  https://reactome.org/download/current/interactors/
#                  reactome.homo_sapiens.interactions.tab-delimited.txt  (human)
#                  reactome.all_species.interactions.tab-delimited.txt   (all organisms)
#                Its pairs are derived by Reactome from complexes and reactions (two proteins in
#                the same complex or the same reaction); they are not curated pairwise interactions,
#                and for non-human organisms the reactions are projected from human by orthology.
#   --gene-list  needed when the *_ids columns hold gene symbols (human): either
#                  https://rest.kegg.jp/list/hsa   or the HGNC "complete set" TSV (hgnc_complete_set.txt).
#                Fetched from rest.kegg.jp when omitted.
#   --out FILE   write the updated table here instead of in place.
#   --date DATE  date to record (default: today).
#
# Each run replaces the earlier entry of the same resource in `evidence`; entries are separated by ' | ':
#   KEGG hsa04110 2026-11-17: relation PPrel:inhibition/phosphorylation
#   KEGG hsa04110 2026-11-17: both on map, no relation
#   KEGG hsa04110 2026-11-17: not on map: HDAC1
#   Reactome 2026-11-17: complex (R-HSA-69191, R-HSA-69195 +3) PMID 3
#   Reactome 2026-11-17: no pair between source and target
require 'date'
require 'net/http'
require 'optparse'
require 'rexml/document'
require 'set'
require 'uri'

KEGG = 'https://rest.kegg.jp'

def fetch(url, retries: 3)
  retries.times do |i|
    res = Net::HTTP.get_response(URI(url))
    return res.body.force_encoding('UTF-8') if res.is_a?(Net::HTTPSuccess)
    raise "HTTP #{res.code} for #{url}" if i == retries - 1
  rescue StandardError
    raise if i == retries - 1
    sleep 2 * (i + 1)
  end
end

def read_tsv(path)
  lines = File.readlines(path, chomp: true, encoding: 'UTF-8').reject { |l| l.strip.empty? }
  fields = lines.shift.split("\t", -1)
  [fields, lines.map { |l| fields.zip(l.split("\t", -1)).to_h }]
end

def write_tsv(path, fields, rows)
  File.open(path, 'w', encoding: 'UTF-8') do |f|
    f.puts fields.join("\t")
    rows.each { |r| f.puts(r.is_a?(Hash) ? fields.map { |k| r[k] }.join("\t") : r.join("\t")) }
  end
end

# ---- gene symbols -> KEGG gene ids (Entrez Gene ids for human) ----
def symbol_map(text)
  head = text.each_line.first.to_s.chomp.split("\t")
  if head.include?('hgnc_id')   # HGNC complete set
    s, e = head.index('symbol'), head.index('entrez_id')
    return text.each_line.drop(1).to_h { |l| f = l.chomp.split("\t", -1); [f[s], f[e]] }.reject { |_, v| v.to_s.empty? }
  end
  primary = {}
  alias_of = {}
  text.each_line(chomp: true) do |line|   # rest.kegg.jp/list/<org>
    f = line.split("\t")
    id = f[0].to_s.split(':', 2)[1]
    names = f[-1].to_s.split(';', 2)[0].to_s.split(',').map(&:strip).reject(&:empty?)
    next if id.nil? || names.empty?

    primary[names[0]] ||= id
    names.drop(1).each { |n| alias_of[n] ||= id }
  end
  alias_of.merge(primary)
end

# ---- KEGG ----
# Returns [pathway id, entries {entry id => Set[gene id]}, relations, all gene ids]
def parse_kgml(xml_text, org)
  root = REXML::Document.new(xml_text).root
  entries = {}
  root.elements.each('entry') do |e|
    ids = e.attributes['name'].to_s.split.filter_map { |n| n.split(':', 2)[1] if n.start_with?("#{org}:") }
    entries[e.attributes['id']] = ids.to_set
  end
  root.elements.each('entry') do |e|
    next unless e.attributes['type'] == 'group'

    e.elements.each('component') { |c| entries[e.attributes['id']].merge(entries.fetch(c.attributes['id'], [])) }
  end
  rels = root.elements.to_a('relation').map do |r|
    [r.attributes['entry1'], r.attributes['entry2'], r.attributes['type'],
     r.elements.to_a('subtype').map { |s| s.attributes['name'] }]
  end
  [root.attributes['name'].to_s.sub('path:', ''), entries, rels, entries.values.reduce(Set.new, :|)]
end

def kegg_pairs(entries, rels)
  pairs = Hash.new { |h, k| h[k] = Set.new }
  rels.each do |e1, e2, type, subs|
    label = "#{type}:#{subs.empty? ? '-' : subs.join('/')}"
    entries.fetch(e1, []).each do |a|
      entries.fetch(e2, []).each { |b| pairs[Set[a, b]] << label unless a == b }
    end
  end
  pairs
end

def kegg_verdict(src, tgt, name_of, entries, pairs, genes)
  if entries.each_value.any? { |ids| src.any? { |g| ids.include?(g) } && tgt.any? { |g| ids.include?(g) } }
    return 'same entry or complex'
  end

  hits = src.product(tgt).reduce(Set.new) { |acc, (a, b)| pairs.key?(Set[a, b]) ? acc | pairs[Set[a, b]] : acc }
  return "relation #{hits.sort.join('; ')}" unless hits.empty?

  missing = (src + tgt).uniq.reject { |g| genes.include?(g) }.map { |g| name_of.fetch(g, g) }.sort
  missing.empty? ? 'both on map, no relation' : "not on map: #{missing.join(', ')}"
end

# ---- Reactome ----
# Stream the interaction file once and keep the pairs whose two sides both carry an id we ask about.
# Columns: protein id 1, Ensembl gene id(s) 1, Entrez Gene id(s) 1, protein id 2, Ensembl 2, Entrez 2,
#          interaction type, interaction context, PubMed references.
def reactome_pairs(path, wanted)
  bare = ->(field) { field.to_s.split('|').map { |v| v.split(':').last.to_s.strip }.reject(&:empty?) }
  pairs = Hash.new { |h, k| h[k] = { types: Set.new, contexts: Set.new, pmids: Set.new } }
  File.foreach(path, chomp: true, encoding: 'UTF-8') do |line|
    next if line.start_with?('#') || line.empty?

    f = line.split("\t", -1)
    next if f.size < 7

    a = (bare.(f[0]) + bare.(f[1]) + bare.(f[2])).select { |v| wanted.include?(v) }
    next if a.empty?

    b = (bare.(f[3]) + bare.(f[4]) + bare.(f[5])).select { |v| wanted.include?(v) }
    next if b.empty?

    a.product(b).each do |x, y|
      next if x == y

      rec = pairs[Set[x, y]]
      rec[:types] << f[6].to_s.strip unless f[6].to_s.strip.empty?
      rec[:contexts].merge(bare.(f[7]))
      rec[:pmids].merge(bare.(f[8]).grep(/\A\d+\z/))
    end
  end
  pairs
end

def reactome_verdict(src, tgt, pairs)
  hit = src.product(tgt).filter_map { |a, b| pairs[Set[a, b]] if pairs.key?(Set[a, b]) }
  return 'no pair between source and target' if hit.empty?

  types = hit.flat_map { |h| h[:types].to_a }.uniq.sort
  ctx = hit.flat_map { |h| h[:contexts].to_a }.uniq.sort
  pmids = hit.flat_map { |h| h[:pmids].to_a }.uniq
  shown = ctx.first(3).join(', ') + (ctx.size > 3 ? " +#{ctx.size - 3}" : '')
  "#{types.join('/')} (#{shown}) PMID #{pmids.size}"
end

def put_evidence(row, key, text)
  kept = row['evidence'].to_s.split(' | ').map(&:strip).reject { |e| e.empty? || e.start_with?("#{key} ") }
  row['evidence'] = (kept + [text]).sort.join(' | ')
end

def main
  opt = { kgml: [], date: Date.today.iso8601 }
  parser = OptionParser.new do |o|
    o.banner = 'Usage: check.rb cellcycle_<org>.tsv [--kgml FILE]... [--reactome FILE] [options]'
    o.on('--kgml FILE', 'KGML file of a pathway of this organism (repeatable)') { |v| opt[:kgml] << v }
    o.on('--reactome FILE', 'Reactome interactions, tab-delimited') { |v| opt[:reactome] = v }
    o.on('--gene-list FILE', 'rest.kegg.jp/list/<org> or HGNC complete set, for symbol lookup') { |v| opt[:gene_list] = v }
    o.on('--out FILE', 'write here instead of updating the table in place') { |v| opt[:out] = v }
    o.on('--date DATE', 'date to record (default: today)') { |v| opt[:date] = Date.iso8601(v).iso8601 }
  end
  parser.parse!
  path = ARGV[0]
  abort parser.help if path.nil? || (opt[:kgml].empty? && opt[:reactome].nil?)
  org = File.basename(path)[/\Acellcycle_([a-z]{3,4})\.tsv\z/, 1] or abort 'the table must be named cellcycle_<org>.tsv'

  fields, rows = read_tsv(path)
  %w[source_ids target_ids].each { |c| abort "#{path}: no #{c} column" unless fields.include?(c) }
  fields << 'evidence' unless fields.include?('evidence')
  maps = opt[:kgml].map { |f| parse_kgml(File.read(f, encoding: 'UTF-8'), org) }
  maps.each { |pw, *| abort "#{pw} is not a pathway of #{org}" unless pw.start_with?(org) }

  # translate the ids of the table when they are symbols rather than KEGG gene ids
  ids = rows.flat_map { |r| r['source_ids'].to_s.split('|') + r['target_ids'].to_s.split('|') }.uniq
  on_maps = maps.map(&:last).reduce(Set.new, :|)
  sym = {}
  if org != 'sce' && ids.none? { |g| on_maps.include?(g) }
    sym = symbol_map(opt[:gene_list] ? File.read(opt[:gene_list], encoding: 'UTF-8') : fetch("#{KEGG}/list/#{org}"))
    warn "#{ids.count { |g| sym.key?(g) }}/#{ids.size} gene symbols translated to gene ids"
    lost = ids.reject { |g| sym.key?(g) }
    warn "not translated: #{lost.join(' ')}" unless lost.empty?
  end
  name_of = ids.to_h { |g| [sym.fetch(g, g), g] }
  side = ->(r, col) { r[col].to_s.split('|').map { |g| sym.fetch(g, g) } }

  maps.each do |pw, entries, rels, genes|
    pairs = kegg_pairs(entries, rels)
    tally = Hash.new(0)
    rows.each do |r|
      verdict = kegg_verdict(side.(r, 'source_ids'), side.(r, 'target_ids'), name_of, entries, pairs, genes)
      put_evidence(r, "KEGG #{pw}", "KEGG #{pw} #{opt[:date]}: #{verdict}")
      tally[verdict[/\A(relation|same entry or complex|both on map, no relation|not on map)/, 1]] += 1
    end
    warn "KEGG #{pw} (#{genes.size} genes, #{rels.size} relations): #{tally.map { |k, v| "#{k} #{v}" }.join(', ')}"
  end

  if opt[:reactome]
    wanted = (ids + ids.map { |g| sym.fetch(g, g) }).to_set
    pairs = reactome_pairs(opt[:reactome], wanted)
    found = 0
    rows.each do |r|
      src = (side.(r, 'source_ids') + r['source_ids'].to_s.split('|')).uniq
      tgt = (side.(r, 'target_ids') + r['target_ids'].to_s.split('|')).uniq
      verdict = reactome_verdict(src, tgt, pairs)
      found += 1 unless verdict.start_with?('no pair')
      put_evidence(r, 'Reactome', "Reactome #{opt[:date]}: #{verdict}")
    end
    warn "Reactome: pair found for #{found}/#{rows.size} rows"
  end

  write_tsv(opt[:out] || path, fields, rows)
  warn "updated #{opt[:out] || path}"
end

main if $PROGRAM_NAME == __FILE__
