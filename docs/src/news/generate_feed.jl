# ---------------------------------------------------------------
# RSS feed generator for LowLevelFEM.jl News section
# ---------------------------------------------------------------
using Dates

newsdir = joinpath(@__DIR__)
feedfile = joinpath(newsdir, "feed.xml")

items = String[]

for file in sort(readdir(newsdir))
    if endswith(file, ".md") && file != "index.md"
        path = joinpath(newsdir, file)
        lines = readlines(path)

        # --- Extract title ---
        title = first(filter(x -> startswith(x, "# "), lines))
        title = replace(title, r"^#\s*" => "")

        # --- Extract first 10 non-empty lines as summary ---
        bodylines = filter(x -> !isempty(strip(x)) && !startswith(x, "#"), lines)
        shortbody = join(first(bodylines, min(10, length(bodylines))), " ")

        # --- Escape HTML characters ---
        function esc(s)
            s = replace(s, "&" => "&amp;")
            s = replace(s, "<" => "&lt;")
            s = replace(s, ">" => "&gt;")
            s
        end

        shortbody = esc(shortbody)

        # --- File date & URL ---
        date = first(split(file, "_"))  # e.g. 2025-11-12
        url  = "https://perebalazs.github.io/LowLevelFEM.jl/dev/news/$(file[1:end-3]).html"

        push!(items, """
        <item>
          <title>$(title)</title>
          <link>$(url)</link>
          <description><![CDATA[
          <p>$(shortbody)</p>
          <p><a href="$(url)">Read more...</a></p>
          ]]></description>
          <pubDate>$(date)</pubDate>
        </item>
        """)
    end
end

# --- Assemble full RSS ---
rss = """
<?xml version="1.0" encoding="UTF-8" ?>
<rss version="2.0">
  <channel>
    <title>LowLevelFEM.jl – News</title>
    <link>https://perebalazs.github.io/LowLevelFEM.jl/news/</link>
    <description>Announcements, examples, and updates for the LowLevelFEM Julia package</description>
    <language>en</language>
    <lastBuildDate>$(Dates.format(now(), dateformat"yyyy-mm-dd"))</lastBuildDate>
    $(join(items, "\\n"))
  </channel>
</rss>
"""

write(feedfile, rss)
println("✅ RSS feed generated at: $feedfile")

