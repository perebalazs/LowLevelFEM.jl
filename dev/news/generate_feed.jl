# docs/src/news/generate_feed.jl
using Dates

newsdir = joinpath(@__DIR__)
feedfile = joinpath(newsdir, "feed.xml")

items = String[]

for file in sort(readdir(newsdir))
    if endswith(file, ".md") && file != "index.md"
        path = joinpath(newsdir, file)
        lines = readlines(path)
        title = first(filter(x -> startswith(x, "# "), lines))
        title = replace(title, r"^#\s*" => "")
        date = first(split(file, "_"))  # e.g. 2025-11-12
        url = "https://perebalazs.github.io/LowLevelFEM.jl/news/$(file[1:end-3]).html"
        push!(items, """
        <item>
          <title>$(title)</title>
          <link>$(url)</link>
          <pubDate>$(date)</pubDate>
        </item>
        """)
    end
end

rss = """
<?xml version="1.0" encoding="UTF-8" ?>
<rss version="2.0">
  <channel>
    <title>LowLevelFEM.jl – News</title>
    <link>https://perebalazs.github.io/LowLevelFEM.jl/news/</link>
    <description>Announcements, examples, and updates for the LowLevelFEM Julia package</description>
    <language>en</language>
    <lastBuildDate>$(Dates.format(now(), dateformat"yyyy-mm-dd"))</lastBuildDate>
    $(join(items, "\n"))
  </channel>
</rss>
"""

write(feedfile, rss)
println("✅ RSS feed generated at: $feedfile")

