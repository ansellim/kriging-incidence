-- Run citeproc here (the YAML sets citeproc: false so this filter sees its output),
-- then replace en and em dashes, such as those in page ranges, with hyphens.
local function fix(el)
  local text = el.text:gsub("\u{2013}", "-"):gsub("\u{2014}", "-")
  if text ~= el.text then return pandoc.Str(text) end
end

function Pandoc(doc)
  if doc.meta.bibliography then doc = pandoc.utils.citeproc(doc) end
  return doc:walk({ Str = fix })
end
