#!/usr/bin/env ruby
# frozen_string_literal: true

# Build cellcycle.json for index.html from the master tables in this directory.
# Run it after every edit of a TSV:
#
#   ruby build.rb
#
# Master tables (all tab-separated, UTF-8, first line = header):
#   cellcycle_<org>.tsv  one interaction per row; <org> is a KEGG organism code (hsa, sce, ...).
#                        Human (hsa) is the reference; every other organism is shown against it.
#   orthologs.tsv        hsa_gene, org, gene, relation (ortholog | functional), evidence, curation
#   diseases.tsv         hsa_gene, note_ja, note_en, evidence, curation
#
# Columns of cellcycle_<org>.tsv:
#   id                 stable identifier, unique across all files
#   module             free-text grouping
#   phase_from/to      one of: G1, G1/S, S, G2/M, M-meta, M-ana, M-exit
#   source_genes/ids   acting side: symbols and database ids, '|'-separated, same count
#   target_genes/ids   receiving side, same layout
#   interaction_type   free text; the words inhibition / activation / phosphorylation ... pick the symbol
#   description_ja/en  one or two sentences
#   references         entries separated by ' | '. Each entry is a citation label followed by one optional tag:
#                        <pmid:12345>   PMID confirmed
#                        <https://...>  link to the paper (DOI, PMC ...)
#                        <checked>      citation details confirmed, no identifier recorded
#                      no tag = citation still to be confirmed
#   shown_in           organism(s) in which the cited experiments were done, KEGG codes, '|'-separated
#                      (e.g. xla|hsa: shown in Xenopus egg extract and in human cells). This can differ from
#                      the organism of the table: a human row may rest on frog or fission-yeast experiments.
#   evidence           written by check.rb: what was checked against which resource and when
#   curation           TODO or DONE: manual confirmation of the row against the literature
require 'json'
require 'date'

DIR    = __dir__
PHASES = %w[G1 G1/S S G2/M M-meta M-ana M-exit].freeze
IX_COLUMNS = %w[id module phase_from phase_to source_genes source_ids target_genes target_ids interaction_type
                description_ja description_en references shown_in evidence curation].freeze

def read_tsv(path, required)
  lines = File.readlines(path, chomp: true, encoding: 'UTF-8').reject { |l| l.strip.empty? }
  fields = lines.shift.split("\t", -1)
  missing = required - fields
  abort "#{File.basename(path)}: missing column(s) #{missing.join(', ')}" unless missing.empty?
  lines.each_with_index.map do |l, i|
    row = fields.zip(l.split("\t", -1)).to_h.transform_values { |v| v.to_s.strip }
    row['_at'] = "#{File.basename(path)}:#{i + 2}"
    row
  end
end

def classify(type)
  case type
  when /inhibition|repression|ubiquitination|proteolysis/ then 'inh'
  when /activation|expression|release|recruitment|stabilization/ then 'act'
  when /phosph|regulation/ then 'mod'
  else 'bind'
  end
end

def parse_refs(text)
  text.split(' | ').map(&:strip).reject(&:empty?).map do |entry|
    label, tag = entry.match(/\A(.*?)(?:\s*<([^<>]+)>)?\z/).captures
    ref = { label: label.strip }
    case tag
    when /\Apmid:(\d+)\z/i then ref.merge!(pmid: Regexp.last_match(1), checked: true)
    when %r{\Ahttps?://} then ref.merge!(url: tag, checked: true)
    when 'checked' then ref[:checked] = true
    when nil then ref[:checked] = false
    else raise "unknown reference tag <#{tag}>"
    end
    ref
  end
end

errors = []
seen = {}
orgs = {}
files = Dir[File.join(DIR, 'cellcycle_*.tsv')].sort_by { |f| f.end_with?('_hsa.tsv') ? '' : f }
abort 'cellcycle_hsa.tsv not found' unless files.first&.end_with?('_hsa.tsv')

files.each do |path|
  org = File.basename(path, '.tsv').sub('cellcycle_', '')
  genes = {}
  ix = read_tsv(path, IX_COLUMNS).filter_map do |r|
    at = r['_at']
    src, sid, tgt, tid = %w[source_genes source_ids target_genes target_ids].map { |k| r[k].split('|') }
    from, to = PHASES.index(r['phase_from']), PHASES.index(r['phase_to'])
    problems = []
    problems << 'empty id' if r['id'].empty?
    problems << "duplicate id (also at #{seen[r['id']]})" if seen[r['id']]
    problems << 'unknown phase' unless from && to
    problems << 'phase_to is earlier than phase_from' if from && to && to < from
    problems << 'gene and id counts differ' unless src.size == sid.size && tgt.size == tid.size
    problems << 'no source or target gene' if src.empty? || tgt.empty?
    shown = r['shown_in'].split('|')
    problems << 'shown_in must be KEGG organism codes separated by |' unless shown.any? && shown.all? { |c| c.match?(/\A[a-z]{3,4}\z/) }
    problems << 'curation must be TODO or DONE' unless %w[TODO DONE].include?(r['curation'])
    refs = begin
      parse_refs(r['references'])
    rescue StandardError => e
      problems << e.message
      []
    end
    unless problems.empty?
      errors << "#{at} #{r['id']}: #{problems.join('; ')}"
      next
    end
    seen[r['id']] = at
    src.zip(sid).each { |g, i| genes[g] ||= { id: i } }
    tgt.zip(tid).each { |g, i| genes[g] ||= { id: i } }
    { id: r['id'], module: r['module'], src: src, tgt: tgt, type: r['interaction_type'], cls: classify(r['interaction_type']),
      desc: { ja: r['description_ja'], en: r['description_en'] }, refs: refs, shown: shown, ph: [from, to],
      evidence: r['evidence'].split(' | ').map(&:strip).reject(&:empty?), curation: r['curation'] }
  end
  orgs[org] = { genes: genes, ix: ix }
end

map = read_tsv(File.join(DIR, 'orthologs.tsv'), %w[hsa_gene org gene relation]).filter_map do |r|
  unless %w[ortholog functional].include?(r['relation'])
    errors << "#{r['_at']}: relation must be ortholog or functional"
    next
  end
  { hsa: r['hsa_gene'], org: r['org'], gene: r['gene'], rel: r['relation'] }
end
disease = read_tsv(File.join(DIR, 'diseases.tsv'), %w[hsa_gene note_ja note_en])
          .to_h { |r| [r['hsa_gene'], { ja: r['note_ja'], en: r['note_en'] }] }

abort "Not built. Fix these rows first:\n  #{errors.join("\n  ")}" unless errors.empty?

data = { built: Date.today.iso8601, phases: PHASES, orgs: orgs, map: map, disease: disease }
File.write(File.join(DIR, 'cellcycle.json'), JSON.generate(data), encoding: 'UTF-8')

orgs.each do |org, o|
  done = o[:ix].count { |x| x[:curation] == 'DONE' }
  refs = o[:ix].flat_map { |x| x[:refs] }
  warn format('%<org>s: %<n>d interactions, %<g>d genes, curation DONE %<d>d/%<n>d, references confirmed %<c>d/%<r>d, rows with check.rb evidence %<e>d',
              org: org, n: o[:ix].size, g: o[:genes].size, d: done, c: refs.count { |q| q[:checked] }, r: refs.size,
              e: o[:ix].count { |x| x[:evidence].any? })
end
unknown = map.reject { |m| orgs[m[:org]]&.dig(:genes, m[:gene]) }.map { |m| "#{m[:org]}:#{m[:gene]}" }
warn "orthologs.tsv: #{map.size} pairs#{unknown.empty? ? '' : "; not used in any interaction: #{unknown.join(' ')}"}"
warn 'wrote cellcycle.json'
